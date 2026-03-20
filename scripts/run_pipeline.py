#!/usr/bin/env python3
"""
run_pipeline.py
---------------
Command-line entry point for the AS TWAS + COLOC pipeline.

Usage
-----
  # Demo / mock run (no external tools required)
  python scripts/run_pipeline.py --mock

  # Real run (requires GTEx data, MetaXcan, R + coloc)
  python scripts/run_pipeline.py --config config/as.yaml

Options
-------
  --config PATH   Path to YAML config file (default: config/as.yaml)
  --mock          Run in mock/demo mode (default: True if --no-mock not set)
  --no-mock       Run with real external tools (requires full setup)
  --base-dir DIR  Repository root directory (default: current directory)
"""

import argparse
import sys
from pathlib import Path

# Add repository root to sys.path so ``src`` is importable
_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_REPO_ROOT))

from src.pipeline import run_pipeline
from src.logging_utils import get_logger

logger = get_logger("run_pipeline")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="AS TWAS + COLOC Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config",
        default=None,
        help="Path to YAML configuration file (default: config/as.yaml)",
    )
    mock_group = parser.add_mutually_exclusive_group()
    mock_group.add_argument(
        "--mock",
        action="store_true",
        default=True,
        help="Run in mock/demo mode (default)",
    )
    mock_group.add_argument(
        "--no-mock",
        dest="mock",
        action="store_false",
        help="Run with real external tools (requires full data setup)",
    )
    parser.add_argument(
        "--base-dir",
        default=".",
        help="Repository root directory (default: current directory)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    logger.info(
        "Starting pipeline (mock=%s, config=%s, base_dir=%s)",
        args.mock,
        args.config,
        args.base_dir,
    )

    try:
        summary = run_pipeline(
            config_path=args.config,
            mock=args.mock,
            base_dir=args.base_dir,
        )
        logger.info("Pipeline finished successfully.")
        logger.info("Summary: %s", summary)
        return 0
    except Exception as exc:
        logger.error("Pipeline failed: %s", exc, exc_info=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
