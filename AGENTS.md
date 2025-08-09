# AGENTS.md — Developer Agents & CI Automation

This repo uses lightweight agents (scripts + LLM prompts) to automate routine tasks: code review hints, documentation updates, mesh QA, and benchmark triage. Agents run locally via CLI and in CI on pull requests.

> These are *assistive* and never auto‑commit to `main` without human review. They operate under strict permission boundaries and produce artifacts for maintainers to inspect.

## Agent Catalog

### 1) Mesh QA Agent
**Purpose:** Validate Gmsh meshes (or imported STEP→mesh) for quality: min dihedral angles, aspect ratios, boundary tags, and port cross‑section sanity.
- **Inputs:** `.msh`/STEP, mesh config, thresholds
- **Outputs:** HTML/JSON report; annotated screenshots; suggested fixes
- **Triggers:** `tools/agents/mesh_qa.py --in examples/*.msh`

### 2) Numerics QA Agent
**Purpose:** Spot numerical risks in PRs (ill‑conditioned matrices, missing complex arithmetic, curl‑curl sign conventions, boundary application).
- **Inputs:** Diff context, unit test logs
- **Outputs:** Markdown review comments with citations to code lines
- **Triggers:** CI step on PR; local pre‑commit hook

### 3) Docs Agent
**Purpose:** Keep `README.md`, API docs, and tutorials consistent with code changes.
- **Inputs:** Public headers, CLI help, example scripts
- **Outputs:** PR comment with proposed doc patches
- **Triggers:** `tools/agents/docs_sync.py` or CI job

### 4) Benchmark Agent
**Purpose:** Run validation pack nightly and on PRs touching core; compare against tolerances.
- **Inputs:** Benchmark YAML, reference results
- **Outputs:** Trend plots, pass/fail summary; regression artifacts
- **Triggers:** `tools/bench/run_pack.py --suite golden`

### 5) Release Notes Agent
**Purpose:** Generate `CHANGELOG` candidate from merged PRs and labels.
- **Inputs:** Git log, PR labels
- **Outputs:** Markdown section for the upcoming version
- **Triggers:** `tools/agents/release_notes.py`

## Implementation
Agents are simple Python CLIs. Optionally they can call an LLM via your provider of choice (e.g., Bedrock, OpenAI) or a **local MCP server** for context retrieval.

````

tools/
agents/
mesh\_qa.py
numerics\_review\.py
docs\_sync.py
release\_notes.py
bench/
run\_pack.py

````

## CI Wiring (GitHub Actions)
Example CI job that runs build, tests, mesh QA, and benchmark smoke:

```yaml
name: CI
on: [push, pull_request]
jobs:
  build-test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with: { python-version: '3.11' }
      - name: Install deps
        run: |
          sudo apt-get update && sudo apt-get install -y cmake ninja-build
          python -m pip install -r tools/requirements.txt
      - name: Configure
        run: cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
      - name: Build
        run: cmake --build build -j
      - name: Unit tests
        run: ctest --test-dir build -j || (cat build/Testing/Temporary/LastTest.log; exit 1)
      - name: Mesh QA (agents)
        run: python tools/agents/mesh_qa.py --in examples --out artifacts/mesh_qa
      - name: Benchmark Smoke
        run: python tools/bench/run_pack.py --suite smoke --out artifacts/bench
      - name: Upload artifacts
        uses: actions/upload-artifact@v4
        with:
          name: vectorem-artifacts
          path: artifacts
````

### Optional macOS local smoke

On Apple Silicon (M3 Max):

```bash
tools/scripts/setup_macos.sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build -L smoke -j
```

## Optional: MCP Server & RAG KB

If you use an MCP server for repository context and retrieval:

```json
{
  "mcp": {
    "servers": {
      "repo-kb": {
        "command": "python",
        "args": ["tools/mcp/repo_kb.py"],
        "env": {"KB_PATH": "./docs"}
      }
    }
  }
}
```

Then agents call MCP tools to fetch docs/snippets for grounded suggestions.

## Security & Guardrails

* Agents run with read‑only tokens; write changes as patches in PR comments
* Large file access and network calls are restricted in CI
* Logs and artifacts are retained for audit; secrets are masked by Actions

## Local Usage Examples

```bash
python tools/agents/mesh_qa.py --in examples/cube_cavity.msh
python tools/agents/numerics_review.py --diff HEAD~1..HEAD
python tools/bench/run_pack.py --suite golden --filter patch_2p45GHz
```

## Adding a New Agent

1. Create `tools/agents/<name>.py` with `--in`, `--out`,args and JSON output
2. Add a Make/CI step and README snippet
3. Provide a sample input and expected output in `tools/agents/samples/`

```
