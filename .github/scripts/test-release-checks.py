"""Exercise the release gate without contacting GitHub or dispatching jobs."""

import json
import os
from pathlib import Path
import subprocess
import tempfile


workflow = Path('.github/workflows/release.yaml').read_text()
step = workflow.split('      - name: Require successful main checks\n', 1)[1]
block = step.split('        run: |\n', 1)[1].split('\n      - ', 1)[0]
body = '\n'.join(line[10:] if line.startswith('          ') else line
                 for line in block.splitlines())
required = [
    'ubuntu-latest (release)', 'windows-latest (release)', 'macos-latest (release)',
    'generated-files', 'metadata', 'coverage', 'minimum-testthat',
    'hard-dependencies-only', 'pdf-manual', 'C and C++', 'benchmark',
    'cran-dependency-audit',
]
successful = [dict(name=name, status='completed', conclusion='success',
                   started_at='2026-08-30T01:00:00Z') for name in required]

wrapper = r'''
api_calls=0
gh() {
  case "$1" in
    api)
      api_calls=$((api_calls + 1))
      if [[ "$CASE" == missing || "$CASE" == skipped ]] && (( api_calls > 1 )); then
        cat "$FIXTURE_ROOT/success.json"
      elif [[ "$CASE" == pending ]] && (( api_calls > 2 )); then
        cat "$FIXTURE_ROOT/success.json"
      else
        cat "$FIXTURE_ROOT/initial.json"
      fi
      ;;
    workflow)
      printf '%s\n' "$*" >> "$FIXTURE_ROOT/dispatched"
      ;;
    *) return 90 ;;
  esac
}
sleep() { :; }
source "$FIXTURE_ROOT/body.sh"
'''

for case in ('success', 'missing', 'skipped', 'failed', 'pending', 'old-failure'):
    records = [dict(record) for record in successful]
    if case == 'missing':
        records = [record for record in records if record['name'] != 'windows-latest (release)']
    elif case == 'skipped':
        records[2]['conclusion'] = 'skipped'
    elif case == 'failed':
        records[2]['conclusion'] = 'failure'
    elif case == 'pending':
        records[2].update(status='in_progress', conclusion=None)
    elif case == 'old-failure':
        records.append(dict(name=required[2], status='completed', conclusion='failure',
                            started_at='2026-08-29T01:00:00Z'))
    with tempfile.TemporaryDirectory(prefix='kfl1ou-release-fixture-') as directory:
        root = Path(directory)
        (root / 'body.sh').write_text(body)
        (root / 'initial.json').write_text(json.dumps([dict(check_runs=records)]))
        (root / 'success.json').write_text(json.dumps([dict(check_runs=successful)]))
        env = dict(os.environ, CASE=case, FIXTURE_ROOT=directory, RUNNER_TEMP=directory,
                   GITHUB_REPOSITORY='fixture/project', RELEASE_SHA='a' * 40,
                   RELEASE_TAG='v3.0.6')
        result = subprocess.run(['bash', '-e', '-o', 'pipefail', '-c', wrapper],
                                env=env, capture_output=True, text=True, timeout=15)
        assert (result.returncode == 0) == (case != 'failed'), (case, result.stderr)
        dispatched = (root / 'dispatched').read_text() if (root / 'dispatched').exists() else ''
        expected = 'workflow run R-CMD-check.yaml --ref v3.0.6\n' if case in ('missing', 'skipped') else ''
        assert dispatched == expected, (case, dispatched)

print('Release validation fixtures passed (6 cases).')
