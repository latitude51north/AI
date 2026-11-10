# ABOUTME: Tests that compose.yaml keeps the properties the pattern depends on.
# ABOUTME: Also checks that no token is stored in the repository.

import yaml

from conftest import PROJECT


def services():
    return yaml.safe_load((PROJECT / "compose.yaml").read_text())["services"]


def test_sidecar_shares_the_app_network_namespace():
    proxy = services()["api-proxy"]
    assert proxy["network_mode"] == "service:app"
    assert "app" in proxy["depends_on"]


def test_sidecar_publishes_no_ports_of_its_own():
    # Docker rejects ports on a container that borrows another's network.
    proxy = services()["api-proxy"]
    assert "ports" not in proxy
    assert "networks" not in proxy


def test_sidecar_forwards_to_the_api_on_loopback():
    assert services()["api-proxy"]["command"] == "caddy reverse-proxy --from :8081 --to 127.0.0.1:8080"


def test_proxy_port_is_published_on_one_address_and_defaults_to_loopback():
    assert services()["app"]["ports"] == ["${HOST_ADDRESS:-127.0.0.1}:8081:8081"]


def test_api_port_itself_is_never_published():
    assert not [port for port in services()["app"]["ports"] if "8080" in port]


def test_token_comes_from_the_environment_and_is_required():
    token = services()["app"]["environment"]["API_TOKEN"]
    assert token.startswith("${API_TOKEN:?")


def test_env_file_is_ignored_and_only_a_placeholder_is_committed():
    assert ".env" in (PROJECT / ".gitignore").read_text().splitlines()
    example = (PROJECT / ".env.example").read_text()
    assert "API_TOKEN=replace-with-a-long-random-string" in example
    assert "HOST_ADDRESS=127.0.0.1" in example
