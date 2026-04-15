from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import date, datetime
from pathlib import Path
import json
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen


HAL_API_URL = "https://api.archives-ouvertes.fr/search/"
DEFAULT_HAL_QUERY = "*:*"
DEFAULT_HAL_ROWS = 5
DEFAULT_HAL_SORT = "producedDate_tdate desc"
DEFAULT_HAL_COLLECTIONS = ("FEEL", "CEMOSIS")
DEFAULT_HAL_FIELDS = (
    "docid",
    "title_s",
    "uri_s",
    "producedDate_tdate",
    "producedDateY_i",
    "docType_s",
    "journalTitle_s",
    "bookTitle_s",
    "authFullName_s",
    "doiId_s",
)

HAL_DOC_TYPE_LABELS = {
    "ART": "Journal article",
    "COMM": "Conference paper",
    "OUV": "Book",
    "COUV": "Book chapter",
    "THESE": "Thesis",
    "HDR": "HDR",
    "REPORT": "Report",
    "UNDEFINED": "Preprint",
}


@dataclass(frozen=True)
class HalPublication:
    docid: str
    title: str
    url: str
    produced_at: str | None
    year: int | None
    doc_type: str | None
    venue: str | None
    authors: tuple[str, ...]
    doi: str | None

    def as_dict(self) -> dict[str, object]:
        return asdict(self)

    @property
    def doc_type_label(self) -> str | None:
        if not self.doc_type:
            return None
        return HAL_DOC_TYPE_LABELS.get(self.doc_type, self.doc_type)


def _first_value(value):
    if isinstance(value, list):
        return value[0] if value else None
    return value


def _coerce_since_value(value: str | date | datetime | None) -> str | None:
    if value is None:
        return None
    if isinstance(value, datetime):
        if value.tzinfo is None:
            return value.strftime("%Y-%m-%dT%H:%M:%SZ")
        return value.astimezone().strftime("%Y-%m-%dT%H:%M:%SZ")
    if isinstance(value, date):
        return f"{value.isoformat()}T00:00:00Z"
    text = str(value).strip()
    if not text:
        return None
    if "T" in text:
        return text
    return f"{text}T00:00:00Z"


def _coerce_collections(value) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        parts = [part.strip() for part in value.split(",")]
        return tuple(part for part in parts if part)
    if isinstance(value, (list, tuple)):
        return tuple(str(item).strip() for item in value if str(item).strip())
    return ()


def load_hal_publication_config(repo_root: str | Path) -> dict[str, object]:
    repo_root = Path(repo_root)
    config_path = repo_root / ".github" / "plan-ci.json"
    try:
        payload = json.loads(config_path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return {}
    release_notes = payload.get("releaseNotes", {})
    publications = release_notes.get("publications", {})
    hal = publications.get("hal", {})
    return hal if isinstance(hal, dict) else {}


class HalPublicationService:
    def __init__(self, repo_root: str | Path | None = None) -> None:
        self.repo_root = Path(repo_root or ".").resolve()
        self.config = load_hal_publication_config(self.repo_root)

    def fetch(
        self,
        *,
        query: str | None = None,
        rows: int | None = None,
        sort: str | None = None,
        collections: tuple[str, ...] | list[str] | str | None = None,
        since: str | date | datetime | None = None,
    ) -> list[HalPublication]:
        resolved_collections = (
            _coerce_collections(collections)
            or _coerce_collections(self.config.get("collections"))
            or DEFAULT_HAL_COLLECTIONS
        )
        params = {
            "q": query or str(self.config.get("query") or DEFAULT_HAL_QUERY),
            "rows": str(rows if rows is not None else int(self.config.get("rows") or DEFAULT_HAL_ROWS)),
            "sort": sort or str(self.config.get("sort") or DEFAULT_HAL_SORT),
            "wt": "json",
            "fl": ",".join(DEFAULT_HAL_FIELDS),
        }
        filters = []
        if resolved_collections:
            filters.append(f"collCode_s:({' OR '.join(resolved_collections)})")
        since_value = _coerce_since_value(since or self.config.get("since"))
        if since_value:
            filters.append(f"producedDate_tdate:[{since_value} TO NOW/DAY]")
        if filters:
            params["fq"] = filters
        request = Request(
            f"{HAL_API_URL}?{urlencode(params, doseq=True)}",
            headers={"User-Agent": "feelpp-ops/1.0"},
        )
        try:
            with urlopen(request, timeout=20) as response:
                payload = json.load(response)
        except (HTTPError, URLError) as exc:
            raise RuntimeError(f"HAL lookup failed: {exc}") from exc

        docs = payload.get("response", {}).get("docs", [])
        publications: list[HalPublication] = []
        for doc in docs:
            publications.append(
                HalPublication(
                    docid=str(doc.get("docid") or ""),
                    title=str(_first_value(doc.get("title_s")) or doc.get("label_s") or ""),
                    url=str(doc.get("uri_s") or ""),
                    produced_at=_first_value(doc.get("producedDate_tdate")),
                    year=doc.get("producedDateY_i"),
                    doc_type=_first_value(doc.get("docType_s")),
                    venue=str(_first_value(doc.get("journalTitle_s")) or _first_value(doc.get("bookTitle_s")) or ""),
                    authors=tuple(doc.get("authFullName_s") or ()),
                    doi=str(doc.get("doiId_s") or "") or None,
                )
            )
        return publications

    def format_markdown(
        self,
        publications: list[HalPublication],
        *,
        heading: str = "## Recent Publications using Feel++",
    ) -> str:
        if not publications:
            return ""
        lines = [heading, ""]
        for publication in publications:
            authors = list(publication.authors)
            if len(authors) > 3:
                author_text = ", ".join(authors[:3]) + ", et al."
            else:
                author_text = ", ".join(authors)
            meta = []
            if publication.doc_type_label:
                meta.append(publication.doc_type_label)
            if publication.venue:
                meta.append(publication.venue)
            if publication.produced_at:
                meta.append(str(publication.produced_at)[:10])
            suffix = f" ({'; '.join(meta)})" if meta else ""
            doi_suffix = f", DOI: `{publication.doi}`" if publication.doi else ""
            if author_text:
                lines.append(f"- [{publication.title}]({publication.url}) — {author_text}{suffix}{doi_suffix}")
            else:
                lines.append(f"- [{publication.title}]({publication.url}){suffix}{doi_suffix}")
        return "\n".join(lines)
