#!/usr/bin/env python3
"""
LLM-ready context crawler for the MrBiomics project.

REQUIREMENTS
install required software (tested with gitingest v0.3.1)
```
conda create -n crawler python=3.11 requests pip
pip install gitingest
```

GITHUB TOKEN
Use with personal GitHub token such that GitHub API is not rate limited.
Generate a GitHub personal access token (Classic) from your profile:

  - Go to https://github.com → profile avatar (top-right) → Settings.
  - In the left sidebar choose Developer settings → Personal access tokens → Tokens (classic) → Generate new token (classic).
  - Give it a note (e.g., “MrBiomics crawler”), pick an expiry that fits your needs, and select minimum scopes:
      - repo (full control of private repos) if you need private access.
      - For public-only access you can leave all scopes unchecked; GitHub will still issue a token that raises rate limits.
  - Create token, copy it immediately (GitHub shows it once) and keep it secret.

  Use it by exporting GITHUB_TOKEN="<token>" in your shell before running the script, pass --token <token> on the command line, or add it to your `.bashrc` profile.
"""

# standard libraries
from __future__ import annotations

import argparse
import datetime as dt
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence
from urllib.parse import urlparse

# libraries to be installed
import requests
from gitingest import ingest

MODULE_EXCLUDE_PATTERNS = (
  "test/**",
  "tests/**",
  "docs/**",
  "**/test/**",
  "**/tests/**",
  "**/docs/**",
)


@dataclass
class RepoSpec:
  label: str
  url: str
  include_patterns: Sequence[str] = field(default_factory=tuple)
  exclude_patterns: Sequence[str] = field(default_factory=tuple)

  def __post_init__(self) -> None:
      self.include_patterns = tuple(self.include_patterns)
      self.exclude_patterns = tuple(self.exclude_patterns)


MAIN_REPO = RepoSpec(
  label="MrBiomics Main",
  url="https://github.com/epigen/MrBiomics",
)

MODULE_REPOS = [
  RepoSpec("unsupervised_analysis", "https://github.com/epigen/unsupervised_analysis", exclude_patterns=MODULE_EXCLUDE_PATTERNS),
  RepoSpec("fetch_ngs", "https://github.com/epigen/fetch_ngs", exclude_patterns=MODULE_EXCLUDE_PATTERNS),
  RepoSpec("spilterlize_integrate", "https://github.com/epigen/spilterlize_integrate", exclude_patterns=MODULE_EXCLUDE_PATTERNS),
  RepoSpec("dea_limma", "https://github.com/epigen/dea_limma", exclude_patterns=MODULE_EXCLUDE_PATTERNS),
  RepoSpec("enrichment_analysis", "https://github.com/epigen/enrichment_analysis", exclude_patterns=MODULE_EXCLUDE_PATTERNS),
  RepoSpec("genome_tracks", "https://github.com/epigen/genome_tracks", exclude_patterns=MODULE_EXCLUDE_PATTERNS),
  RepoSpec("atacseq_pipeline", "https://github.com/epigen/atacseq_pipeline", exclude_patterns=MODULE_EXCLUDE_PATTERNS),
  RepoSpec("rnaseq_pipeline", "https://github.com/epigen/rnaseq_pipeline", exclude_patterns=MODULE_EXCLUDE_PATTERNS),
  RepoSpec("scrnaseq_processing_seurat", "https://github.com/epigen/scrnaseq_processing_seurat", exclude_patterns=MODULE_EXCLUDE_PATTERNS),
  RepoSpec("dea_seurat", "https://github.com/epigen/dea_seurat", exclude_patterns=MODULE_EXCLUDE_PATTERNS),
  RepoSpec("mixscape_seurat", "https://github.com/epigen/mixscape_seurat", exclude_patterns=MODULE_EXCLUDE_PATTERNS),
]

WIKI_REPO = RepoSpec(
  label="MrBiomics Wiki",
  url="https://github.com/epigen/MrBiomics.wiki",
)


def repo_identifier(repo_url: str) -> str:
  parsed = urlparse(repo_url)
  path = parsed.path.strip("/")
  if path.endswith(".git"):
      path = path[:-4]
  return path


def github_headers(token: Optional[str]) -> dict[str, str]:
  headers = {
      "Accept": "application/vnd.github+json",
      "User-Agent": "MrBiomics-context-crawler",
  }
  if token:
      headers["Authorization"] = f"Bearer {token}"
  return headers


def fetch_latest_commit(repo_url: str, token: Optional[str] = None) -> tuple[str, str]:
  """Use the GitHub API to look up the default-branch commit for a normal repo."""
  repo_id = repo_identifier(repo_url)
  if not repo_id:
      raise RuntimeError(f"Unable to parse repository identifier from {repo_url}")
  session = requests.Session()
  session.headers.update(github_headers(token))
  try:
      repo_resp = session.get(f"https://api.github.com/repos/{repo_id}", timeout=30)
      repo_resp.raise_for_status()
      default_branch = repo_resp.json().get("default_branch", "main")
      commit_resp = session.get(
          f"https://api.github.com/repos/{repo_id}/commits/{default_branch}",
          timeout=30,
      )
      commit_resp.raise_for_status()
      data = commit_resp.json()
  except requests.HTTPError as exc:
      if exc.response is not None and exc.response.status_code == 404:
          return "unknown", "unknown"
      raise RuntimeError(f"Failed to fetch commit metadata for {repo_url}") from exc
  except requests.RequestException as exc:
      raise RuntimeError(f"Failed to fetch commit metadata for {repo_url}") from exc
  sha = data.get("sha") or "unknown"
  commit_date = data.get("commit", {}).get("committer", {}).get("date") or "unknown"
  return sha, commit_date


def fetch_repo_digest(repo: RepoSpec, token: Optional[str] = None) -> str:
  kwargs: dict[str, object] = {}
  if token:
      kwargs["token"] = token
  if repo.include_patterns:
      kwargs["include_patterns"] = list(repo.include_patterns)
  if repo.exclude_patterns:
      kwargs["exclude_patterns"] = list(repo.exclude_patterns)
  try:
      summary, tree, content = ingest(repo.url, **kwargs)
  except Exception as exc:
      warning = f"GitIngest failed for {repo.url}: {exc}"
      print(f"WARNING: {warning}", file=sys.stderr)
      return warning
  parts = []
  for part in (summary, tree, content):
      if part:
          parts.append(part.strip())
  return "\n\n".join(parts).strip()


def fetch_wiki_via_git(repo: RepoSpec, token: Optional[str] = None) -> tuple[str, str, str]:
  """
  Clone the wiki git repository, concatenate pages with separators, and return
  the digest along with HEAD commit metadata.
  """
  clone_url = repo.url
  if token:
      clone_url = clone_url.replace("https://", f"https://{token}:x-oauth-basic@")

  with tempfile.TemporaryDirectory() as tmpdir:
      try:
          subprocess.run(
              ["git", "clone", "--depth=1", clone_url, tmpdir],
              check=True,
              stdout=subprocess.PIPE,
              stderr=subprocess.PIPE,
          )
      except subprocess.CalledProcessError as exc:
          stderr = exc.stderr.decode("utf-8", "ignore")
          raise RuntimeError(f"Git clone failed for {repo.url}: {stderr}") from exc

      sha = "unknown"
      commit_date = "unknown"
      try:
          result = subprocess.run(
              ["git", "-C", tmpdir, "rev-parse", "HEAD"],
              check=True,
              stdout=subprocess.PIPE,
              stderr=subprocess.PIPE,
              text=True,
          )
          sha = result.stdout.strip() or "unknown"
      except subprocess.CalledProcessError:
          pass

      try:
          result = subprocess.run(
              ["git", "-C", tmpdir, "log", "-1", "--format=%cI"],
              check=True,
              stdout=subprocess.PIPE,
              stderr=subprocess.PIPE,
              text=True,
          )
          commit_date = result.stdout.strip() or "unknown"
      except subprocess.CalledProcessError:
          pass

      pages: list[str] = []
      root = Path(tmpdir)
      for path in sorted(root.rglob("*")):
          if not path.is_file():
              continue
          rel = path.relative_to(root)
          if rel.name == ".gitignore":
              continue
          if ".git" in rel.parts:
              continue
          text = path.read_text(encoding="utf-8", errors="replace").strip()
          sections = [
              "================================================================",
              f"WIKI PAGE: {rel}",
              "================================================================",
              "",
              text,
          ]
          pages.append("\n".join(sections).strip())

  digest = "\n\n".join(pages) if pages else "No wiki pages found."
  return digest, sha, commit_date


def build_single_repo_file(repo: RepoSpec, generated_at: dt.datetime, token: Optional[str]) -> str:
  sha, commit_date = fetch_latest_commit(repo.url, token)
  digest = fetch_repo_digest(repo, token)
  lines = [
      f"Snapshot generated: {generated_at.isoformat()}",
      f"Repository label: {repo.label}",
      f"Repository URL: {repo.url}",
      f"Latest commit SHA: {sha}",
      f"Latest commit date: {commit_date}",
      "Captured via GitIngest (summary, directory tree, file contents).",
      "",
      digest,
  ]
  return "\n".join(lines).rstrip() + "\n"


def build_modules_file(modules: Sequence[RepoSpec], generated_at: dt.datetime, token: Optional[str]) -> str:
  header_lines = [
      f"Snapshot generated: {generated_at.isoformat()}",
      "Modules included (latest default branch commits):",
  ]
  for repo in modules:
      header_lines.append(f"  - {repo.label} ({repo.url})")
  header_lines.extend(
      [
          "",
          'Each module section begins with "===== MODULE: <label> =====".',
      ]
  )
  sections = []
  for repo in modules:
      sha, commit_date = fetch_latest_commit(repo.url, token)
      digest = fetch_repo_digest(repo, token)
      section_lines = [
          "================================================================",
          f"===== MODULE: {repo.label} =====",
          f"Repository URL: {repo.url}",
          f"Latest commit SHA: {sha}",
          f"Latest commit date: {commit_date}",
          "Excluded patterns: " + (", ".join(repo.exclude_patterns) if repo.exclude_patterns else "None"),
          "================================================================",
          "",
          digest,
      ]
      sections.append("\n".join(section_lines).strip())
  header_text = "\n".join(header_lines).strip()
  return "\n\n".join([header_text] + sections).rstrip() + "\n"


def build_wiki_file(repo: RepoSpec, generated_at: dt.datetime, token: Optional[str]) -> str:
  digest, sha, commit_date = fetch_wiki_via_git(repo, token)
  lines = [
      f"Snapshot generated: {generated_at.isoformat()}",
      f"Repository label: {repo.label}",
      f"Repository URL: {repo.url}",
      f"Latest commit SHA: {sha}",
      f"Latest commit date: {commit_date}",
      "Captured via git clone of the wiki repository; each page is delimited below.",
      "",
      digest,
  ]
  return "\n".join(lines).rstrip() + "\n"


def prepare_output_dir(root: Path, dry_run: bool) -> Path:
  folder_name = f"MrBiomics_AI_resources"
  target_dir = (root / folder_name).resolve()
  if dry_run:
      print(f"[dry-run] Would create directory {target_dir}")
  else:
      target_dir.mkdir(parents=True, exist_ok=True)
  return target_dir


def write_text(path: Path, content: str, dry_run: bool) -> None:
  if dry_run:
      print(f"[dry-run] Would write {path} ({len(content)} characters)")
      return
  path.write_text(content, encoding="utf-8")
  print(f"Wrote {path} ({len(content)} characters)")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
  parser = argparse.ArgumentParser(description="Collect LLM-ready MrBiomics resources via GitIngest.")
  parser.add_argument("--output-root", type=Path, default=Path.cwd(), help="Directory where the resources folder will be created.")
  parser.add_argument("--token", type=str, default=None, help="GitHub token for higher rate limits when calling the GitHub API and GitIngest.")
  parser.add_argument("--skip-main", action="store_true", help="Skip the main MrBiomics repository snapshot.")
  parser.add_argument("--skip-modules", action="store_true", help="Skip the module repositories snapshot.")
  parser.add_argument("--skip-wiki", action="store_true", help="Skip the MrBiomics wiki snapshot.")
  parser.add_argument("--dry-run", action="store_true", help="Print the planned actions without writing files.")
  return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
  if ingest is None:
      raise SystemExit("GitIngest is required; install it with `pip install gitingest` before running this script.")
  args = parse_args(argv)
  generated_at = dt.datetime.now(dt.timezone.utc)
  output_root = args.output_root.resolve()
  output_dir = prepare_output_dir(output_root, args.dry_run)
  try:
      if not args.skip_main:
          main_text = build_single_repo_file(MAIN_REPO, generated_at, args.token)
          write_text(output_dir / "MrBiomics_main.txt", main_text, args.dry_run)
      if not args.skip_modules:
          modules_text = build_modules_file(MODULE_REPOS, generated_at, args.token)
          write_text(output_dir / "MrBiomics_modules.txt", modules_text, args.dry_run)
      if not args.skip_wiki:
          wiki_text = build_wiki_file(WIKI_REPO, generated_at, args.token)
          write_text(output_dir / "MrBiomics_wiki.txt", wiki_text, args.dry_run)
  except RuntimeError as exc:
      print(f"ERROR: {exc}", file=sys.stderr)
      return 1
  if args.dry_run:
      print(f"[dry-run] Completed without writing files. Target directory: {output_dir}")
  else:
      print(f"Completed. Resources available in {output_dir}")
  return 0


if __name__ == "__main__":
  sys.exit(main())

