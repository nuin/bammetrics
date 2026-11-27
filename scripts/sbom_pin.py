#!/usr/bin/env python3
"""Make a CycloneDX SBOM reproducible: set metadata.timestamp and serialNumber.

    sbom_pin.py SBOM.json COMMIT_SHA SOURCE_DATE_EPOCH

The timestamp is SOURCE_DATE_EPOCH in RFC 3339 UTC. The serial number is a
UUID v5 of the commit's URL, so one commit always gives one SBOM.
"""
import json
import sys
import uuid
from datetime import datetime, timezone

COMMIT_URL = "https://github.com/nuin/bammetrics/commit/"


def pin(text, sha, epoch):
    doc = json.loads(text)
    stamp = datetime.fromtimestamp(int(epoch), timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    doc["serialNumber"] = "urn:uuid:" + str(uuid.uuid5(uuid.NAMESPACE_URL, COMMIT_URL + sha))
    doc["metadata"]["timestamp"] = stamp
    return json.dumps(doc, indent=2, ensure_ascii=False)


def main(argv):
    if len(argv) != 4:
        print("usage: sbom_pin.py SBOM.json COMMIT_SHA SOURCE_DATE_EPOCH", file=sys.stderr)
        return 2
    path, sha, epoch = argv[1:]
    with open(path, encoding="utf-8") as f:
        text = f.read()
    with open(path, "w", encoding="utf-8") as f:
        f.write(pin(text, sha, epoch))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
