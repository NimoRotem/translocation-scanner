"""Event aggregation and throttling for SSE streaming.

Sits between the pipeline and SSE subscribers. Batches high-frequency events
and rate-limits to maintain ≤50 events/sec to each client.
"""
from __future__ import annotations
import asyncio
import time
import logging
from collections import defaultdict
from typing import Optional

logger = logging.getLogger(__name__)


class EventAggregator:
    """Aggregates pipeline events into throttled SSE-friendly batches.

    - pair.density: batched every 500ms
    - chrom.bin_update: batched every 250ms
    - evidence.highlight: rate-limited to 20/sec, biased split > clip > discordant
    - Total cap: 50 events/sec
    """

    def __init__(self, max_events_per_sec: int = 50):
        self.max_events_per_sec = max_events_per_sec
        self._subscribers: dict[str, asyncio.Queue] = {}
        self._density_buffer: dict[tuple, dict] = {}
        self._bin_buffer: dict[tuple, dict] = {}
        self._highlight_buffer: list[dict] = []
        self._passthrough_buffer: list[dict] = []
        self._last_density_flush = time.monotonic()
        self._last_bin_flush = time.monotonic()
        self._last_highlight_flush = time.monotonic()
        self._events_this_second = 0
        self._second_start = time.monotonic()
        self._lock = None  # initialized lazily in async context
        self._event_counter = 0

    def _get_lock(self):
        if self._lock is None:
            self._lock = asyncio.Lock()
        return self._lock

    def subscribe(self, client_id: str) -> asyncio.Queue:
        """Register a new SSE client."""
        q: asyncio.Queue = asyncio.Queue(maxsize=500)
        self._subscribers[client_id] = q
        logger.info("SSE client subscribed: %s (total: %d)", client_id, len(self._subscribers))
        return q

    def unsubscribe(self, client_id: str):
        """Remove an SSE client."""
        self._subscribers.pop(client_id, None)
        logger.info("SSE client unsubscribed: %s (total: %d)", client_id, len(self._subscribers))

    def push_event(self, event: dict):
        """Push a raw pipeline event into the aggregator (called from pipeline thread)."""
        event_type = event.get("type", "")
        self._event_counter += 1
        event["_seq"] = self._event_counter

        if event_type == "pair.density":
            key = (event.get("chrom_a", ""), event.get("chrom_b", ""))
            self._density_buffer[key] = event
        elif event_type == "chrom.bin_update":
            key = (event.get("chrom", ""), event.get("bin_start", 0))
            self._bin_buffer[key] = event
        elif event_type == "evidence.highlight":
            self._highlight_buffer.append(event)
        else:
            # Passthrough: scan.started, scan.completed, scan.stage_changed,
            # scan.progress, scan.throughput, provisional.top_bins,
            # validation.*, chrom.progress, error
            self._passthrough_buffer.append(event)

    async def flush(self):
        """Flush buffered events to all subscribers, respecting rate limits."""
        now = time.monotonic()
        events_to_send = []

        # Reset per-second counter
        if now - self._second_start >= 1.0:
            self._events_this_second = 0
            self._second_start = now

        # Always send passthrough events first (stage changes, progress, etc.)
        while self._passthrough_buffer:
            events_to_send.append(self._passthrough_buffer.pop(0))

        # Flush density buffer every 500ms
        if now - self._last_density_flush >= 0.5 and self._density_buffer:
            for evt in self._density_buffer.values():
                events_to_send.append(evt)
            self._density_buffer.clear()
            self._last_density_flush = now

        # Flush bin updates every 250ms
        if now - self._last_bin_flush >= 0.25 and self._bin_buffer:
            for evt in self._bin_buffer.values():
                events_to_send.append(evt)
            self._bin_buffer.clear()
            self._last_bin_flush = now

        # Flush highlights: rate-limit to 20/sec, prefer split > clip > discordant
        if now - self._last_highlight_flush >= 0.05 and self._highlight_buffer:
            # Sort by priority: split first, then clip, then discordant
            priority = {"split": 0, "clip_pileup": 1, "new_pair_cluster": 2, "discordant": 3}
            self._highlight_buffer.sort(
                key=lambda e: priority.get(e.get("evidence_type", ""), 9)
            )
            # Take up to 20 per second worth
            budget = min(20, self.max_events_per_sec - len(events_to_send))
            if budget > 0:
                events_to_send.extend(self._highlight_buffer[:budget])
                self._highlight_buffer = self._highlight_buffer[budget:]
            # Drop excess if buffer grows too large
            if len(self._highlight_buffer) > 200:
                self._highlight_buffer = self._highlight_buffer[-100:]
            self._last_highlight_flush = now

        # Cap total events
        if len(events_to_send) > self.max_events_per_sec:
            events_to_send = events_to_send[:self.max_events_per_sec]

        # Fan out to subscribers
        dead_clients = []
        for client_id, queue in self._subscribers.items():
            for evt in events_to_send:
                try:
                    queue.put_nowait(evt)
                except asyncio.QueueFull:
                    # Drop oldest events for slow clients
                    try:
                        queue.get_nowait()
                        queue.put_nowait(evt)
                    except (asyncio.QueueEmpty, asyncio.QueueFull):
                        pass

        for cid in dead_clients:
            self.unsubscribe(cid)

        self._events_this_second += len(events_to_send)
        return len(events_to_send)

    async def run_flush_loop(self, interval: float = 0.1):
        """Run the flush loop — call this as an asyncio task."""
        while True:
            try:
                await self.flush()
            except Exception:
                logger.exception("Error in aggregator flush loop")
            await asyncio.sleep(interval)

    def has_subscribers(self) -> bool:
        return len(self._subscribers) > 0
