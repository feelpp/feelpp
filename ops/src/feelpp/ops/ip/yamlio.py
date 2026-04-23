from __future__ import annotations

import json
from pathlib import Path
from typing import Any

try:
    import yaml
except ModuleNotFoundError:  # pragma: no cover - exercised in minimal test envs
    yaml = None


def _line_indent(line: str) -> int:
    return len(line) - len(line.lstrip(" "))


def _strip_continuations(lines: list[str]) -> list[str]:
    normalized: list[str] = []
    for line in lines:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        stripped = line.strip()
        if normalized and ":" not in stripped and stripped != "-" and not stripped.startswith("- "):
            normalized[-1] = f"{normalized[-1]} {stripped}"
            continue
        normalized.append(line.rstrip())
    return normalized


def _parse_scalar(raw: str) -> Any:
    value = raw.strip()
    if value in {"null", "~"}:
        return None
    if value == "[]":
        return []
    if value == "{}":
        return {}
    if value == "true":
        return True
    if value == "false":
        return False
    if value.isdigit():
        return int(value)
    if (value.startswith('"') and value.endswith('"')) or (
        value.startswith("'") and value.endswith("'")
    ):
        return value[1:-1]
    return value


def _next_content_indent(lines: list[str], index: int) -> int | None:
    if index >= len(lines):
        return None
    return _line_indent(lines[index])


def _parse_simple_yaml_block(lines: list[str], index: int, indent: int) -> tuple[Any, int]:
    if index >= len(lines):
        return None, index
    current = lines[index]
    current_indent = _line_indent(current)
    if current_indent < indent:
        return None, index
    stripped = current.strip()
    if stripped == "-" or stripped.startswith("- "):
        values: list[Any] = []
        while index < len(lines):
            line = lines[index]
            line_indent = _line_indent(line)
            text = line.strip()
            if line_indent != indent or (text != "-" and not text.startswith("- ")):
                break
            item_text = "" if text == "-" else text[2:].strip()
            index += 1
            if not item_text:
                child_indent = _next_content_indent(lines, index)
                if child_indent is None:
                    values.append(None)
                    continue
                child, index = _parse_simple_yaml_block(lines, index, child_indent)
                values.append(child)
                continue
            if ":" in item_text:
                key, raw_value = item_text.split(":", 1)
                item: dict[str, Any] = {}
                if raw_value.strip():
                    item[key.strip()] = _parse_scalar(raw_value)
                else:
                    child_indent = _next_content_indent(lines, index)
                    if child_indent is None:
                        item[key.strip()] = None
                    else:
                        child, index = _parse_simple_yaml_block(lines, index, child_indent)
                        item[key.strip()] = child
                if index < len(lines) and _line_indent(lines[index]) > indent:
                    child, index = _parse_simple_yaml_block(lines, index, _line_indent(lines[index]))
                    if isinstance(child, dict):
                        item.update(child)
                    else:
                        raise ValueError("Expected YAML mapping after list mapping item")
                values.append(item)
                continue
            values.append(_parse_scalar(item_text))
        return values, index

    mapping: dict[str, Any] = {}
    while index < len(lines):
        line = lines[index]
        line_indent = _line_indent(line)
        if line_indent < indent:
            break
        if line_indent != indent:
            raise ValueError(f"Unexpected YAML indentation: {line}")
        text = line.strip()
        if text.startswith("- "):
            break
        key, separator, raw_value = text.partition(":")
        if not separator:
            raise ValueError(f"Expected YAML mapping entry: {line}")
        index += 1
        if raw_value.strip():
            mapping[key.strip()] = _parse_scalar(raw_value)
            continue
        child_indent = _next_content_indent(lines, index)
        if child_indent is None:
            mapping[key.strip()] = None
            continue
        if child_indent < indent:
            mapping[key.strip()] = None
            continue
        child, index = _parse_simple_yaml_block(lines, index, child_indent)
        mapping[key.strip()] = child
    return mapping, index


def _load_simple_yaml(text: str) -> dict[str, Any]:
    try:
        payload = json.loads(text)
    except json.JSONDecodeError:
        lines = _strip_continuations(text.splitlines())
        if not lines:
            return {}
        payload, index = _parse_simple_yaml_block(lines, 0, _line_indent(lines[0]))
        if index != len(lines):
            raise ValueError("Unable to parse complete YAML document")
    if not isinstance(payload, dict):
        raise ValueError("Expected YAML mapping")
    return payload


def _render_scalar(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    return json.dumps(str(value), ensure_ascii=False)


def _dump_simple_yaml(value: Any, indent: int = 0) -> list[str]:
    prefix = " " * indent
    if isinstance(value, dict):
        if not value:
            return [f"{prefix}{{}}"]
        lines: list[str] = []
        for key, child in value.items():
            if child == []:
                lines.append(f"{prefix}{key}: []")
                continue
            elif child == {}:
                lines.append(f"{prefix}{key}: {{}}")
                continue
            if isinstance(child, (dict, list)):
                lines.append(f"{prefix}{key}:")
                lines.extend(_dump_simple_yaml(child, indent + 2))
            else:
                lines.append(f"{prefix}{key}: {_render_scalar(child)}")
        return lines
    if isinstance(value, list):
        if not value:
            return [f"{prefix}[]"]
        lines = []
        for item in value:
            if isinstance(item, (dict, list)):
                lines.append(f"{prefix}-")
                lines.extend(_dump_simple_yaml(item, indent + 2))
            else:
                lines.append(f"{prefix}- {_render_scalar(item)}")
        return lines
    return [f"{prefix}{_render_scalar(value)}"]


def load_yaml(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    payload = yaml.safe_load(text) if yaml is not None else _load_simple_yaml(text)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected YAML mapping in {path}")
    return payload


def dump_yaml(payload: dict[str, Any]) -> str:
    if yaml is None:
        return "\n".join(_dump_simple_yaml(payload)) + "\n"
    return yaml.safe_dump(
        payload,
        allow_unicode=True,
        default_flow_style=False,
        sort_keys=False,
    )


def write_yaml(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(dump_yaml(payload), encoding="utf-8")
