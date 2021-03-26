import http.server
import socketserver
import webbrowser
import json
import shutil
import logging

from functools import partial
from collections import defaultdict

import numpy as np
import scipy
from scipy.cluster.hierarchy import linkage

from cblaster.classes import Session
from cblaster.helpers import get_project_root
from cblaster.extract_clusters import get_sorted_cluster_hierarchies


LOG = logging.getLogger(__name__)


def transform_linkage_matrix(matrix):
    """Converts SciPy linkage matrix to D3 hierarchical format."""

    hierarchy = {}
    total = matrix.shape[0] + 1  # Linkage matrix is n-1 by 4

    # Generate placeholders to pop for each label
    for index in range(total):
        hierarchy[index] = {"name": index}

    for index, (one, two, distance, count) in enumerate(matrix):
        one = int(one)
        two = int(two)
        new = total + index  # total rows + index
        hierarchy[new] = {
            "name": new,
            "length": distance,
            "children": [hierarchy.pop(one), hierarchy.pop(two)]
        }

    return hierarchy[new]


def generate_linkage_matrix(array):
    """Generate a normalised linkage matrix from a given array."""
    array = [
        [cell["value"] for cell in cells]
        for cells in array
    ]
    array = np.array(array)
    array = array / np.max(array)
    return linkage(array, "ward")


def get_cell(query, cluster, cluster_id):
    hits = [
        {
            "name": hit.subject,
            "identity": hit.identity,
            "coverage": hit.coverage,
            "bitscore": hit.bitscore,
            "evalue": hit.evalue,
            "strand": subject.strand,
            "start": subject.start,
            "end": subject.end,
            "ipg": subject.ipg,
        }
        for subject in cluster
        for hit in subject.hits
        if hit.query == query
    ]
    if not hits:
        value = 0
    elif not hits[0].get("identity", None):
        value = len(hits) * 100
    else:
        value = max(hit["identity"] for hit in hits)
    cell = {
        "query": query,
        "cluster": cluster_id,
        "value": value,
        "hits": hits,
        "flag": -1,
    }
    return cell


def flag_duplicate_cells(cells):
    """Finds cells that contain hits present in other cells."""
    groups = defaultdict(list)
    for index, cell in enumerate(cells):
        for hit in cell["hits"]:
            name = hit["name"]
            groups[name].append(index)

    # Assign numbers to each unique shared hit group.
    # Sorts groups, then checks for index overlap. If none, iterate number.
    group, number = set(), 0
    for indices in sorted(groups.values(), key=min):
        if len(indices) <= 1:
            continue
        if not group or not group.isdisjoint(indices):
            group.update(indices)
        else:
            number += 1
        for index in indices:
            cells[index]["flag"] = number


def get_data(session, sort_clusters=False, max_clusters=None):
    matrix = []
    labels = {}
    counts = {
        "queries": len(session.queries),
        "hits": 0,
        "subjects": 0,
        "clusters": 0,
        "scaffolds": 0,
        "organisms": 0
    }
    if sort_clusters:
        cluster_hierarchy = get_sorted_cluster_hierarchies(
            session,
            max_clusters=max_clusters
        )
    else:
        cluster_hierarchy = [
            (cluster, scaffold, organism.full_name)
            for organism in session.organisms
            for accession, scaffold in organism.scaffolds.items()
            for cluster in scaffold.clusters
        ]
    cluster_id = 0
    organisms = set()
    scaffolds = set()
    for cluster, scaffold, organism_name in cluster_hierarchy:
        organisms.add(organism_name)
        scaffolds.add(scaffold.accession)
        counts["subjects"] += len(scaffold.subjects)
        counts["hits"] += sum(len(sub.hits) for sub in scaffold.subjects)
        counts["clusters"] += 1

        # Save the cluster name and scaffold
        labels[cluster_id] = {
            "id": cluster_id,
            "name": organism_name,
            "scaffold": scaffold.accession,
            "start": cluster.start,
            "end": cluster.end,
            "number": cluster.number,
            "score": round(cluster.score, 3)
        }

        # Generate all cells for the heatmap
        cells = [
            get_cell(query, cluster, cluster_id)
            for query in session.queries
        ]

        # Flag cells which contain hits present in other cells
        flag_duplicate_cells(cells)
        matrix.append(cells)
        cluster_id += 1

    counts["organisms"] = len(organisms)
    counts["scaffolds"] = len(scaffolds)

    # Only generate a linkage matrix if there is more than one result
    if len(matrix) > 1 and not sort_clusters:
        linkage_matrix = generate_linkage_matrix(matrix)
        hierarchy = transform_linkage_matrix(linkage_matrix)
    elif sort_clusters:
        inner = [{"name": number, "children": []} for number in range(0, counts["clusters"])]
        hierarchy = {"name": counts["clusters"], "children": inner}
    else:
        hierarchy = {
            "name": 0,
            "length": 0,
            "children": []
        }

    return {
        "queries": session.queries,
        "labels": labels,
        "counts": counts,
        "matrix": matrix,
        "hierarchy": hierarchy,
        "sort_clusters": sort_clusters
    }


class CustomHandler(http.server.BaseHTTPRequestHandler):
    """Handler for serving cblaster plots."""

    def __init__(self, data, chart, *args, **kwargs):
        self._data = data
        self._chart = chart
        self._dir = get_project_root() / "plot"
        super().__init__(*args, **kwargs)

    def copy_file(self, source):
        shutil.copyfileobj(source, self.wfile)

    def send_headers(self, mime):
        self.send_response(200)
        self.send_header("Content-Type", mime)
        self.end_headers()

    def log_message(self, format, *args):
        """Suppresses logging messages on every request."""
        return

    def do_GET(self):
        """Serves each component of the cblaster plot."""
        if self.path == "/data.json":
            self.send_headers("text/json")
            self.wfile.write(json.dumps(self._data).encode())
            return
        path, mime = None, None
        if self.path == "/":
            if self._chart == "heatmap":
                path, mime = self._dir / "cblaster.html", "text/html"
            elif self._chart == "gne":
                path, mime = self._dir / "gne.html", "text/html"
        elif self.path == "/index.css":
            path, mime = self._dir / "index.css", "text/css"
        elif self.path == "/d3.min.js":
            path, mime = self._dir / "d3.min.js", "text/javascript"
        elif self.path == "/cblaster.js":
            path, mime = self._dir / "cblaster.js", "text/javascript"
        elif self.path == "/gne.js":
            path, mime = self._dir / "gne.js", "text/javascript"
        if not path:
            return
        with path.open("rb") as fp:
            self.send_headers(mime)
            self.copy_file(fp)


def save_html(data, output, chart="heatmap"):
    """Generates a static HTML file with all visualisation code."""

    if chart == "heatmap":
        base, script = "cblaster.html", "cblaster.js"
    elif chart == "gne":
        base, script = "gne.html", "gne.js"
    else:
        raise ValueError("Invalid chart specified, expected 'heatmap' or 'gne'")

    directory = get_project_root() / "plot"

    with (directory / base).open() as fp:
        html = fp.read()

    css_string = '<link href="index.css" rel="stylesheet"></link>'
    d3_string = '<script src="d3.min.js"></script>'
    cb_string = f'<script src="{script}"></script>'

    with (directory / "index.css").open() as fp:
        css = fp.read()
        html = html.replace(css_string, f"<style>{css}</style>")

    with (directory / "d3.min.js").open() as fp:
        d3 = fp.read()
        html = html.replace(d3_string, f"<script>{d3}</script>")

    with (directory / script).open() as fp:
        js = f"const data={json.dumps(data)}" + fp.read()
        html = html.replace(cb_string, f"<script>{js}</script>")

    with open(output, "w") as fp:
        fp.write(html)


def serve_html(data, chart="heatmap"):
    handler = partial(CustomHandler, data, chart)

    # Instantiate a new server, bind to any open port
    with socketserver.TCPServer(("localhost", 0), handler) as httpd:

        # Automatically open web browser to bound address
        address, port = httpd.server_address
        url = f"http://{address}:{port}/"
        webbrowser.open(url)

        # Start serving the plot; shutdown on a keyboard interrupt
        try:
            LOG.info(f"Serving cblaster plot at {url} (Ctrl+C to stop).")
            httpd.serve_forever()
        except KeyboardInterrupt:
            httpd.shutdown()


def plot_session(session, output=None, sort_clusters=False, max_clusters=None, testing=False):
    data = get_data(session, sort_clusters, max_clusters=max_clusters)
    if output:
        LOG.info(f"Saving cblaster plot HTML to: {output}")
        save_html(data, output)
        if not testing:
            webbrowser.open(output)
    else:
        serve_html(data)


def plot_gne(data, output=None, testing=False):
    if output:
        LOG.info(f"Saving gne plot HTML to: {output}")
        save_html(data, chart="gne", output=output)
        if not testing:
            webbrowser.open(output)
    else:
        serve_html(data, chart="gne")


def plot_session_file(path, output=None):
    with open(path) as fp:
        session = Session.from_json(fp)
    plot_session(session, output=output)
