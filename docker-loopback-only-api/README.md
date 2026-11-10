# Reaching a loopback-only API inside a Docker container

The code for the post
[Reaching a loopback-only API inside a Docker container](https://latitude51north.com/blog/docker-loopback-only-api/).

Some services listen only on `127.0.0.1` and have no setting to change that. Inside a container
this puts them out of reach. A proxy on the host sends to the host's own loopback. A published
port delivers to the container's network address, where nothing is listening.

A second container that shares the first one's network namespace can reach the service, because
its `127.0.0.1` is the service's `127.0.0.1`. The post explains why. This folder is a version you
can run.

Every address, port and name here is an example. Addresses starting `192.0.2` belong to the
range reserved for documentation.

## What is here

| File | What it does |
| --- | --- |
| `compose.yaml` | The app, and a [Caddy](https://caddyserver.com/docs/command-line#caddy-reverse-proxy) sidecar in the app's network namespace. |
| `app/server.py` | A stand-in for the app. Its API listens on `127.0.0.1:8080` only and requires a bearer token. |
| `firewall.sh` | Limits the published port to one source address, in Docker's `DOCKER-USER` chain. |
| `verify.sh` | The three checks from the post. |
| `.env.example` | The two settings, with placeholder values. |
| `tests/` | Tests for all of the above. |

The stand-in is there so the example runs by itself. To use the pattern for real, keep the
`api-proxy` service and the `ports` entry, and put your own app where `app` is. Change `8080` in
the sidecar's command to the port your service listens on.

## Run it

You need Docker with Compose.

```bash
cp .env.example .env
```

Set `API_TOKEN` in `.env` to a long random string. `openssl rand -hex 32` prints one. Then start
both containers:

```bash
docker compose up -d
```

As shipped, the proxy's port is published on the host's loopback, so only the host can reach it.
Check it from the host:

```bash
TOKEN=$(grep '^API_TOKEN=' .env | cut -d= -f2-) ./verify.sh http://127.0.0.1:8081
```

```text
PASS  without the token: 401
PASS  with the token: 200
```

The 401 is the more useful line. It proves the whole path works, from the client through the
published port and the sidecar to the API, and that the API still refuses a caller with no token.

## Let one other machine in

The proxy makes a localhost-only API reachable from the network. Three controls replace the
protection the loopback address gave it.

**Publish on one address.** In `.env`, set `HOST_ADDRESS` to the host's own network address, such
as `192.0.2.10`, and run `docker compose up -d` again. The port is then published on that address
and no other.

**Allow one source.** Published ports bypass ufw, so the rule goes in the `DOCKER-USER` chain. Print
the rules first, then apply them:

```bash
./firewall.sh --dry-run 192.0.2.20
sudo ./firewall.sh 192.0.2.20
```

```text
iptables -I DOCKER-USER -p tcp --dport 8081 -j DROP
iptables -I DOCKER-USER -p tcp --dport 8081 -s 192.0.2.20 -j ACCEPT
iptables -L DOCKER-USER -n --line-numbers
```

`192.0.2.20` is the one machine allowed in. The script ends by listing the chain. ACCEPT must be
above DROP. It is safe to run again: it does not add a second DROP, and it puts ACCEPT back on top.

**Keep the token.** The proxy does not add the token. Caddy passes the caller's `Authorization`
header through and the API checks it. A caller needs both an allowed address and the token.

Then check from both sides. From the allowed machine:

```bash
TOKEN=<the token> ./verify.sh http://192.0.2.10:8081
```

From any other machine, where the request should get no answer:

```bash
./verify.sh --blocked http://192.0.2.10:8081
```

## What to know before relying on it

- **The firewall rules match every container.** `DOCKER-USER` sees a request after its destination
  has been rewritten, so the rules match the port inside the container, for all forwarded traffic.
  Choose a port nothing else uses. `firewall.sh` takes the port as a second argument.
- **The rules do not survive a reboot.** On Debian and Ubuntu, install `iptables-persistent` and
  run `sudo netfilter-persistent save`.
- **Changing the allowed address leaves the old rule behind.** `firewall.sh` adds an ACCEPT for the
  address you give it. It does not remove one added earlier for a different address. Delete that
  with `sudo iptables -D DOCKER-USER -p tcp --dport 8081 -s <old address> -j ACCEPT`.
- **Docker's nftables mode has no `DOCKER-USER` chain.** There you add the same restriction in a
  table of your own.
- **The sidecar depends on the app's container.** It borrows that container's network. The Compose
  file sets `restart: true` under `depends_on`, so Compose restarts the sidecar when it replaces
  the app's container. That does not cover a container that dies and is restarted by Docker. If
  the proxy stops answering, run `docker compose restart api-proxy`.
- **It serves plain HTTP.** The token crosses the network unencrypted between the client and the
  proxy. Give `--from` a hostname and Caddy serves HTTPS. `--internal-certs` makes it issue the
  certificate from its own certificate authority, which the client then has to trust.

## Tests

```bash
uv run --with pytest --with pyyaml pytest tests/ -q
```

The tests start the stand-in API as a real process and call it over HTTP, including through
`verify.sh`. They check that `compose.yaml` keeps the properties the pattern depends on.
`firewall.sh` is run against a recording stand-in for `iptables`.

The tests do not start Docker, and they never change a real firewall rule.
