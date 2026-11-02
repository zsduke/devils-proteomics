"""English abbreviated-month display; stored timestamps stay ISO."""
from datetime import datetime
import re

_TIMESTAMP = re.compile(r'^(\d{4})-(\d{2})-(\d{2})(?:[T ](\d{2}:\d{2}:\d{2}(?:\.\d+)?)(Z|[+-]\d{2}:?\d{2})?)?$')
_TIMESTAMP_FIELDS = frozenset({'created_utc', 'created_at', 'downloaded_utc', 'generated_utc', 'retrieved_at', 'attempted_at'})
_MONTHS = ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')


def display_timestamp(value):
    """Reorder a valid ISO date without converting its recorded clock or offset."""
    if not isinstance(value, str):
        return value
    match = _TIMESTAMP.fullmatch(value)
    if match is None:
        return value
    try:
        datetime.fromisoformat(value)
    except ValueError:
        return value
    year, month, day, clock, offset = match.groups()
    return f'{_MONTHS[int(month)-1]} {int(day)}, {year}' + (f' {clock}{offset or ""}' if clock else '')


def display_provenance(value):
    """Copy provenance for HTML display, leaving unrelated strings untouched."""
    if isinstance(value, dict):
        return {key: display_timestamp(item) if key in _TIMESTAMP_FIELDS else display_provenance(item)
                for key, item in value.items()}
    if isinstance(value, list):
        return [display_provenance(item) for item in value]
    return value
