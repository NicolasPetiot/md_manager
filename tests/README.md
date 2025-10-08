Commands:

```shell
uv sync              # Will install pytest
uv pip install -e .. # Will install local version od MD-manager

uv run pytest -vv -s    # Perform All tests
uv run pytest -k import # Only perform test_import*.py tests
```
