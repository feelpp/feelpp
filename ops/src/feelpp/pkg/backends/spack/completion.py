# SPDX-FileContributor: Christophe Prud'homme <christophe.prudhomme@feelpp.org>
# SPDX-FileCopyrightText: 2026 University of Strasbourg
# SPDX-License-Identifier: LGPL-3.0-or-later
"""Read-only Bash completion derived from the CLI parser and repository manifests."""
from __future__ import annotations

import argparse
import sys

from feelpp.pkg.cli import build_parser
from .commands.env import _environment_aliases, _environment_manifests, _environment_name
from .runtime import environments_root, resolve_repo_root


def _commands(parser: argparse.ArgumentParser) -> dict[str, argparse.ArgumentParser]:
    """Return child commands without parsing arguments or invoking their handlers."""
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action.choices
    return {}


def _environments(repo_root: str | None) -> set[str]:
    """Offer canonical environment names and aliases that resolve unambiguously."""
    try:
        root = environments_root(resolve_repo_root(repo_root))
        names: set[str] = set()
        aliases: dict[str, set[str]] = {}
        for manifest in _environment_manifests(root):
            name = _environment_name(root, manifest)
            names.add(name)
            for alias in _environment_aliases(root, manifest):
                aliases.setdefault(alias, set()).add(name)
        return names | {alias for alias, owners in aliases.items() if len(owners) == 1}
    except (OSError, ValueError):
        return set()


def complete(words: list[str]) -> list[str]:
    """Complete argument words, including the possibly empty current word.

    Follow the parser tree without validating partial arguments. Remember option
    values so flags do not consume positional environment names. The shell's
    enable command delegates to activate, whose options it also accepts.
    """
    parser = _commands(build_parser())["spack"]
    current = words[-1] if words else ""
    pending = None
    values: dict[str, str] = {}
    path: list[str] = []
    positional_used = False
    for word in words[:-1]:
        if pending is not None:
            values[pending.dest] = word
            pending = None
            continue
        if word.startswith("-"):
            flag, separator, value = word.partition("=")
            action = parser._option_string_actions.get(flag)
            if action is not None and action.nargs != 0:
                if separator:
                    values[action.dest] = value
                else:
                    pending = action
            continue
        command = "activate" if not path and word == "enable" else word
        child = _commands(parser).get(command)
        if child is not None:
            parser = child
            path.append(command)
            positional_used = False
        else:
            positional_used = True

    if pending is not None:
        candidates = (_environments(values.get("repo_root")) if pending.dest == "env"
                      else {str(value) for value in pending.choices} if pending.choices else set())
    elif current.startswith("-"):
        candidates = set(parser._option_string_actions)
    elif _commands(parser):
        candidates = set(_commands(parser))
    elif not positional_used and path in (["activate"], ["env", "show"]):
        candidates = _environments(values.get("repo_root"))
    else:
        candidates = set(parser._option_string_actions)
    return sorted(value for value in candidates if value.startswith(current))


if __name__ == "__main__":
    for candidate in complete(sys.argv[1:]):
        print(candidate)
