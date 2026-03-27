from pathlib import Path

_TOOLBOXES_ROOT = Path(__file__).resolve().parents[2]
_TOOLBOX_ALIASES = {
    "cfpdes": "coefficientformpdes",
}


def toolbox_case(casefile: str) -> str:
    case_path = Path(casefile)
    if not case_path.parts:
        return str(case_path)
    toolbox_name = _TOOLBOX_ALIASES.get(case_path.parts[0], case_path.parts[0])
    return str(_TOOLBOXES_ROOT / toolbox_name / "cases" / Path(*case_path.parts[1:]))
