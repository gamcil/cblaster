from cblaster import parsers


def test_get_arguments_remote_defaults(monkeypatch):
    def mock_path(file_path, access):
        return "test"

    monkeypatch.setattr(parsers, "full_path", mock_path)

    parsed = parsers.parse_args(["search", "-qf", "test"])
    expected = {
        "subcommand": "search",
        "output": None,
        "output_hide_headers": False,
        "output_delimiter": None,
        "output_decimals": 4,
        "binary": None,
        "binary_hide_headers": False,
        "binary_delimiter": None,
        "binary_decimals": 4,
        "debug": False,
        "query_ids": None,
        "query_file": "test",
        "mode": "remote",
        "session_file": None,
        "database": ["nr"],
        "entrez_query": None,
        "rid": None,
        "gap": 20000,
        "unique": 3,
        "min_hits": 3,
        "require": None,
        "max_evalue": 0.01,
        "min_identity": 30,
        "min_coverage": 50,
    }

    for key, value in expected.items():
        assert getattr(parsed, key) == value, f"Expected {key}={value}"
