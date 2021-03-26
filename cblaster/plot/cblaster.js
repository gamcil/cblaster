/* cblaster plot
 * Cameron L.M. Gilchrist, 2020.
 */

const constants = {
	"cellHeight": 27,
	"cellWidth": 44,
	"dendroWidth": 140,
	"showScaffolds": true,
	"multiLineLabels": true,
	"showCounts": true,
	"showCountBorders": true,
	"xAxisOnTop": true,
	"sorted": false,
	"fontFamily": 'system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Ubuntu, "Helvetica Neue", Oxygen, Cantarell, sans-serif'
}

if (typeof data === 'undefined') {
	// Load from HTTP server
	const data = d3.json("data.json").then(data => plot(data));
} else {
	plot(data);
}

function serialise(svg) {
	/* Saves the figure to SVG in its current state.
	 * Clones the provided SVG and sets the width/height of the clone to the
	 * bounding box of the original SVG. Thus, downloaded figures will be sized
	 * correctly.
	 * This function returns a new Blob, which can then be downloaded.
	*/
	node = svg.node();
	const xmlns = "http://www.w3.org/2000/xmlns/";
	const xlinkns = "http://www.w3.org/1999/xlink";
	const svgns = "http://www.w3.org/2000/node";
	const bbox = svg.select("g").node().getBBox()

	node = node.cloneNode(true);
	node.setAttribute("width", bbox.width);
	node.setAttribute("height", bbox.height);
	node.setAttributeNS(xmlns, "xmlns", svgns);
	node.setAttributeNS(xmlns, "xmlns:xlink", xlinkns);

	// Adjust x/y of <g> to account for axis/title position.
	// Replaces the transform attribute, so drag/zoom is ignored.
	d3.select(node)
		.select("g")
		.attr("transform", `translate(${Math.abs(bbox.x)}, ${Math.abs(bbox.y)})`)

	const serializer = new window.XMLSerializer;
	const string = serializer.serializeToString(node);
	return new Blob([string], {type: "image/node+xml"});
}

function download(blob, filename) {
	/* Downloads a given blob to filename.
	 * This function appends a new anchor to the document, which points to the
	 * supplied blob. The anchor.click() method is called to trigger the download,
	 * then the anchor is removed.
	*/
	const link = document.createElement("a");
	link.href = URL.createObjectURL(blob);
	link.download = filename;
	document.body.appendChild(link);
	link.click();
	document.body.removeChild(link);
}

function colourBarGradient(colors) {
	/* Creates a lineargradient element used in the colour bar fill.
	 * This function returns the DOM node, which should be appended to a
	 * <defs> node.
	 * Have to use svg: namespace when creating detached elements.
	*/
	const gradient = d3.create("svg:linearGradient")
		.attr("id", "cbar-gradient")
		.attr("x1", "0%")
		.attr("x2", "100%");
	gradient.append("stop")
		.attr("offset", "0%")
		.attr("stop-color", colors(0));
	gradient.append("stop")
		.attr("offset", "100%")
		.attr("stop-color", colors(1));
	return gradient.node();
}

function colourBar() {
	/* Creates a colour bar element and returns its DOM node.
	 * Must have #cbar-gradient defined.
	*/
	const cbar = d3.create("svg:g")
		.attr("transform", "translate(2, 2)");
	cbar.append("rect")
		.style("fill", "url(#cbar-gradient)");
	cbar.append("rect")
		.style("fill", "none")
		.style("stroke", "black")
		.style("stroke-width", "thin");
	cbar.selectAll("rect")
		.attr("width", 150)
		.attr("height", 10);
	cbar.append("text")
		.text("Identity (%)")
		.attr("x", 75)
		.attr("y", 25)
		.attr("text-anchor", "middle")
	cbar.append("text")
		.text("0")
		.attr("x", 0)
		.attr("y", 25)
		.attr("text-anchor", "start")
	cbar.append("text")
		.text("100")
		.attr("x", 150)
		.attr("y", 25)
		.attr("text-anchor", "end")
	cbar.selectAll("text")
		.style("font-family", constants.fontFamily)
		.style("font-size", "12px")
	return cbar.node();
}

function elbow(d) {
	/* Path generator function to draw right angles in the dendrogram. */
	return `
		M${d.source.y},${d.source.x}
		V${d.target.x}
		H${d.target.y}
	`;
}

function getSummaryHTML(counts) {
	/* Generates the HTML content for the count summary.
	 */
	return `<p>Your search of <b>${counts.queries}</b> queries`
	+ ` returned <b>${counts.hits}</b> hits`
	+ ` from <b>${counts.subjects}</b> unique sequences.</p>`
	+ `<p>cblaster detected <b>${counts.clusters}</b> clusters`
	+ ` across <b>${counts.scaffolds}</b> genomic scaffolds`
	+ ` from <b>${counts.organisms}</b> organisms.</p>`;
}

function flattenArray(array) {
	/* Flattens a 2d array to 1d. */
	return array.reduce((flat, next) => flat.concat(next), []);
}

function getTooltipHTML(d, data) {
	/* Generates the HTML content for a cell hovering tooltip.
	 * It provides the name of the query, as well as a table of each hit.
	 */
	let ncbi = "https://www.ncbi.nlm.nih.gov"
	let cluster = data.labels[d.cluster]
	let mapUrl = `${ncbi}/nuccore/${cluster.scaffold}`
		+ `?report=graph&from=${cluster.start}&to=${cluster.end}`
	d.hits.sort((a, b) => b.identity / b.coverage - a.identity / a.coverage)
	return `
	<p class="tooltip-summary">
		<span>${cluster.name}</span>
		<br>
		<a href="${mapUrl}">
			${getScaffoldString(cluster)}
		</a>
	<table class="tooltip-hits">
	<tbody>
		<tr>
			<th><b>Name</b></th>
			<th><b>Start</b></th>
			<th><b>End</b></th>
			<th><b>Identity</b></th>
			<th><b>Coverage</b></th>
			<th><b>Bitscore</b></th>
			<th><b>E-value</b></th>
		</tr>
		${d.hits.map(h => (`
			<tr>
				<td>
					<a href="${ncbi}/protein/${h.name}">${h.name}</a>
				</td>
				<td>${h.strand === "+" ? h.start : h.end}</td>
				<td>${h.strand === "+" ? h.end : h.start}</td>
				<td>${h.identity ? h.identity.toFixed(2) : "-"}</td>
				<td>${h.coverage ? h.coverage.toFixed(2) : "-"}</td>
				<td>${h.bitscore ? h.bitscore.toFixed(2) : "-"}</td>
				<td>${h.evalue ? h.evalue : "-"}</td>
			</tr>
		`)).join("")}
	</tbody>
	</table>
	`
}

function scaleBranchLengths(node, width) {
	/* Scales the branch lengths of the cluster dendrogram.
	 * This requires cluster distances to be normalised (e.g. divide all distances
	 * in SciPy linkage matrix by the maximum distance). A linear scale maps the
	 * normalised distances to the total width of the dendrogram, and then y values
	 * of every link source and target are scaled accordingly. Terminal targets will
	 * have no length attribute and thus default to 0 (i.e. drawn to max width).
	 */
	const y = d3.scaleLinear()
		.domain([1, 0])  // inversed, since thats how Scipy reports distances
		.range([0, width]);
  for (const link of node.links()) {
		link.source.y = y(link.source.data.length)	
		link.target.y = y(link.target.data.length || 0)	
  }
}

function getScaledTree(hierarchy, width, height) {
	// Generate d3 cluster object from data.
	const root = d3.hierarchy(hierarchy);
	const tree = d3.cluster()
		.size([height, width])
		.separation(() => 1)(root);
	scaleBranchLengths(tree, width);
	return tree;
}

function pruneHierarchy(node, labels) {
	/* Prune any children of a nested object whose name is equal to label.
	 * Additionally, if pruning results in an empty children list, remove the
	 * parent.
	*/
	return node.filter(child => {
		if (child.children) {
			child.children = pruneHierarchy(child.children, labels)	
			if (child.children.length === 0 && constants.sorted !== true)
				return false
		}
    return !labels.includes(child.name)
	})
}

function getScaffoldString(data) {
	let scaffold = data.scaffold + ":" + data.start + "-" + data.end
	return constants.cellHeight < 24 ? `    ${scaffold}` : scaffold
}

function plot(data) {
	data.matrix = flattenArray(data.matrix);
	constants.sorted = data.sort_clusters

  const originalData = JSON.parse(JSON.stringify(data))

	// Reset to the original data. Have to make a deep copy here, since update
	// will mutate data
	d3.select("#btn-reset-filters")
		.on("click", () => {
			const copy = JSON.parse(JSON.stringify(originalData))
			update(copy)
		});

	// Base <svg> and <g> elements. All chart elements go inside <g>. Pan/zoom
	// behaviour is attached to the <svg>.
	const plotDiv = d3.select("#plot");
	const svg = plotDiv.append("svg")
		.classed("wrapper-svg", true)
		.attr("id", "root_svg")
		.attr("cursor", "grab")
		.attr("xmlns", "http://www.w3.org/2000/svg");
	const g = svg.append("g").attr("transform", "translate(2,0)");

	// Set up pan/zoom behaviour, and set default pan/zoom position
	const zoom = d3.zoom()
		.scaleExtent([0, 8])
		.on("zoom", () => g.attr("transform", d3.event.transform))
		.on("start", () => svg.attr("cursor", "grabbing"))
		.on("end", () => svg.attr("cursor", "grab"))
	const transform = d3.zoomIdentity
		.translate(20, 50)
		.scale(1.2)
	svg.call(zoom).call(zoom.transform, transform)

	// Set up bandwidth scales for x and y axes, and heatmap color scale
	const x = d3.scaleBand().padding(0.05)
	const y = d3.scaleBand().padding(0.05)
	const colorScale = d3
		.scaleSequential(d3.interpolateBlues)
		.domain([0, 1])

	// Generate colour bar and linear gradient definition
	const def = svg.append("defs");
	def.append(() => colourBarGradient(colorScale));
	svg.append(() => colourBar());

	// cblaster chart skeleton to be populated in update()
	const dendro = g.append("g")
	const heatmap = g.append("g")
	const heatmapX = heatmap.append("g")
		.attr("class", "g-heatmap-xaxis")
	const heatmapY = heatmap.append("g")
		.attr("class", "g-heatmap-yaxis")
	const heatmapBG = heatmap.append("rect")
		.attr("class", "heatmap-bg")
		.attr("transform", "translate(-1, -1)")
		.style("fill", "white")
		.style("stroke", "#31363b")
		.style("stroke-width", "thin");
	const cellGroup = heatmap.append("g")

	// Set up the save SVG button
	d3.select("#btn-save-svg")
		.on("click", () => {
			const blob = serialise(svg)
			download(blob, "cblaster.svg")
		})

	// Populate the search summary
	const summaryHTML = getSummaryHTML(data.counts)
	d3.select("#p-result-summary")
		.html(summaryHTML)

	// Populate tooltip with current cell data, and adjust position to match the
	// cell in the heatmap (ignoring <g> transforms).
	const cellEnter = (d, i, n) => {
		tooltip
			.html(getTooltipHTML(d, data))
		let rect = n[i].getBoundingClientRect()
		let bbox = tooltip.node().getBoundingClientRect()
		let xOffset = rect.width / 2 - bbox.width / 2
		let yOffset = rect.height * 1.2
		tooltip
			.style("left", rect.x + xOffset + "px")
			.style("top", rect.y + yOffset + "px")
		tooltip.transition()
			.duration(100)
			.style("opacity", 1)
			.style("pointer-events", "all")
	}

	// Transition upon entering the tooltip <div>. This cancels out a previously
	// called transition (i.e. delayed transition in tooltipLeave). Also enables
	// pointer events to allow text selection, clicking hyperlinks, etc.
	const tooltipEnter = () => {
		tooltip.transition()
			.duration(0)
			.style("opacity", 1)
			.style("pointer-events", "all")
	}

	// Delayed tooltip transition for either 1) when user has left heatmap cell
	// and does not go into another, or 2) user entered tooltip <div> and has now
	// left it. Hides tooltip after 400ms, and disables pointer events which would
	// swallow pan/zoom events.
	const tooltipLeave = () => {
		tooltip.transition()
			.delay(400)
			.style("opacity", 0)
			.style("pointer-events", "none")
	}

	const tooltip = plotDiv.append("div")
		.classed("tooltip", true)
		.style("opacity", 0)
		.style("pointer-events", "none")
		.style("position", "absolute")
		.style("padding", "5px")
		.on("mouseenter", tooltipEnter)
		.on("mouseleave", tooltipLeave)

	function update(data) {
		let t = d3.transition().duration(400)

		// Update x-axis domain/range based on current query sequences.
		x.domain(data.queries)
			.range([0, data.queries.length * constants.cellWidth])

		// Generate a scaled d3.cluster object based on the current hierarchy.
		const tree = getScaledTree(
			data.hierarchy,
			constants.dendroWidth,
			constants.cellHeight * Object.keys(data.labels).length,
		);

		// Draw 'elbow' <path> elements for each link in the dendrogram.
		dendro.selectAll("path")
			.data(tree.links())
			.join("path")
			.attr("fill", "none")
			.attr("stroke", "#555")
			.attr("stroke-opacity", 0.4)
			.attr("stroke-width", 1.5)
			.attr("d", elbow);

		// Update y-axis domain/range based on the new dendrogram
		const dendroLeaves = tree.leaves()
		y.domain(dendroLeaves.map(l => l.data.name))
			.range([0, dendroLeaves.length * constants.cellHeight])

		// Plot cells in the heatmap.
		// Expects data.matrix to be flattened, then draws each cell given its
		// corresponding query/cluster number.
		const translate = (d) => `translate(${x(d.query)}, ${y(d.cluster)})`

		// Generate the border colour scale for cells that contain shared hits
		let groupMax = d3.max(data.matrix, d => d.flag)
		let borderColors = d3.scaleSequential(d3.interpolateLab("orange", "red"))
			.domain([0, groupMax])

		// Function to update cell background/count text, since this
		// should be called in both the enter and update selections
		const updateCell = cell => {
			cell.selectAll("rect")
				.attr("width", x.bandwidth() - 1)
				.attr("height", y.bandwidth() - 1)
				.attr("fill", d => colorScale(d.value / 100))
				.style("stroke", d => (constants.showCountBorders && d.flag >= 0) ? borderColors(d.flag) : null)
			cell.selectAll("text")
				.text(d => d.hits.length)
				.attr("x", constants.cellWidth / 2)
				.attr("y", constants.cellHeight / 2 + 4)
				.attr("opacity", (constants.showCounts ? 1 : 0))
				.style("font-size", `${0.5 * constants.cellHeight}px`)
				.style("fill", d => d.value < 60 ? colorScale(1) : colorScale(0))
		}

		// Draw heatmap cells.
		// Each cell contains a <rect> for the background, and <text>
		// for the count label. Cells with no hits (count == 0) are not drawn.
		cellGroup.selectAll("g.cell")
			.data(data.matrix, d => `${d.query}-${d.cluster}`)
			.join(
				enter => {
					enter = enter.append("g")
						.filter(d => d.hits.length > 0)
						.attr("transform", translate)
						.attr("class", "cell")
					enter.append("rect")
						.attr("class", "heatmap-cell")
						.style("stroke-width", "1.5px")
					enter.append("text")
						.attr("class", "heatmap-count")
						.style("text-anchor", "middle")
						.style("font-family", constants.fontFamily)
					enter.call(updateCell)
					return enter
				},
				update => update.call(
					update => {
						update.transition(t)
							.attr("transform", translate)
							.call(updateCell)
					}
				)
			)
			.on("mouseenter", cellEnter)
			.on("mouseleave", tooltipLeave)

		// Calculate current heatmap width/height and transition the background.
		const heatWidth = data.queries.length * constants.cellWidth
		const heatHeight = Object.keys(data.labels).length * constants.cellHeight
		heatmapBG.transition(t)
			.attr("width", heatWidth + 1)
			.attr("height", heatHeight + 2)

		// Calculate new d3 axes for x-axis and call transition.
		if (constants.xAxisOnTop) {
			let axisTop = d3.axisTop(x)
				.tickSize(6)
				.tickPadding(6)
			heatmapX.transition(t)
				.call(axisTop)
				.attr("transform", "translate(0, -1)")
		} else {
			let axisBot = d3.axisBottom(x)
				.tickSize(6)
				.tickPadding(6)
			heatmapX.transition(t)
				.call(axisBot)
				.attr("transform", `translate(0, ${heatHeight})`)
		}

		// Rotate query sequence labels and attach click delete behaviour.
		const clickRemoveX = query => {
      // Escape if there would be no queries left after this action
      let newQueries = data.queries.filter(q => q !== query)
      if (newQueries.length === 0)
        return

      // Preliminary filter; remove any cells matching the deleted query label,
      // also mark clusters which still have hits so we can delete empty ones
      // after
      let filteredCells = []
      let filteredClusters = new Set()
      for (const cell of data.matrix) {
        if (cell.query === query) continue
        if (cell.hits.length > 0)
          filteredClusters.add(cell.cluster)
        filteredCells.push(cell)
      }

      // Determine which clusters to keep/remove from the plot
      let newClusters = {}
      let removeClusters = []
      for (let label of Object.keys(data.labels)) {
        label = parseInt(label)
        if (filteredClusters.has(label))
          newClusters[label] = data.labels[label]
        else
          removeClusters.push(label)
      }

      // Escape if there would be no clusters left after this action
      if (Object.keys(newClusters).length === 0)
        return

      // Prune the hierarchy if any labels are removed
      let hierarchy = removeClusters.length > 0
        ? { ...data.hierarchy, children: pruneHierarchy(data.hierarchy.children, removeClusters) }
        : { ...data.hierarchy }

      // Call plot update with the new data
			update({
				...data,
				queries: newQueries,
        matrix: filteredCells.filter(cell => filteredClusters.has(cell.cluster)),
        labels: newClusters,
				hierarchy: hierarchy
      })
		}
		heatmapX.selectAll("text")
			.style("font-family", constants.fontFamily)
			.style("font-size", `${Math.min(constants.cellWidth * 0.35, 14)}px`)
			.style("text-anchor", (constants.xAxisOnTop) ? "start" : "end")
			.attr("dy", (constants.xAxisOnTop) ? "" : ".6em")
			.attr("transform", "rotate(-25)")
			.on("click", clickRemoveX)

		// Convert cluster IDs to multiline text labels with organism name and
		// cluster scaffold coordinates, and attach click delete behaviour.
		const labelScaffold = (g) => {
			g.selectAll("text").remove()
			let text = g.append("text")
			text.append("tspan")
        .text(d => `${data.labels[d].name} (Cluster ${data.labels[d].number}, score: ${data.labels[d].score})`)
				.attr("x", 10)
				.attr("dy", constants.multiLineLabels ? "-0.1em": ".3em")
				.attr("text-anchor", "start")
				.style("fill", "#31363b")
				.style("font-family", constants.fontFamily)
				.style("font-style", "italic")
				.style(
					"font-size",
					constants.multiLineLabels
					? `${Math.min(0.5 * constants.cellHeight, 12)}px`
					: `${Math.min(0.7 * constants.cellHeight, 14)}px`
				)
			text.append("tspan")
				.text(d => getScaffoldString(data.labels[d]))
				.attr("x", constants.multiLineLabels ? 10 : null)
				.attr("dy", constants.multiLineLabels ? "1.1em": null) 
				.attr("text-anchor", "start")
				.style("font-style", "normal")
				.style("font-family", constants.fontFamily)
				.style(
					"font-size",
					constants.multiLineLabels
					? `${Math.min(0.4 * constants.cellHeight, 10)}px`
					: `${Math.min(0.6 * constants.cellHeight, 12)}px`
				)
				.style("fill", "grey")
		}
		const labelNoScaffold = (g) => {
			g.selectAll("text").remove()
			g.append("text")
        .text(d => `${data.labels[d].name} (Cluster ${data.labels[d].number}, score: ${data.labels[d].score})`)
				.attr("x", 10)
				.attr("dy", ".3em")
				.attr("text-anchor", "start")
				.style("fill", "#31363b")
				.style("font-family", constants.fontFamily)
				.style("font-style", "italic")
				.style("font-size", `${0.7 * constants.cellHeight}px`)
		}
		const pickLabel = g => {
			if (constants.showScaffolds) {
				labelScaffold(g)
			} else {
				labelNoScaffold(g)
			}
		}

		// Draw cluster labels and ticks.
		// This is implemented with a separate join to more easily enable
		// multi-line labels and toggles with transitions and click to remove behaviour.
		let yAxisTranslate = d => `translate(${heatWidth}, ${y(d) + y.bandwidth() / 2})`
		let clickRemoveY = label => {
      // Escape if there would be no clusters left after this action
      let newLabels = {}
      for (let lbl of Object.keys(data.labels)) {
        lbl = parseInt(lbl)
        if (lbl !== label) newLabels[lbl] = data.labels[lbl]
      }
      if (Object.keys(newLabels).length === 0)
        return

      // Update plot with new data
			update({
				...data,
				hierarchy: {
					...data.hierarchy,
					"children": pruneHierarchy(data.hierarchy.children, [label])
				},
				matrix: data.matrix.filter(cell => cell.cluster !== label),
        labels: newLabels
      })
		}
		heatmapY.selectAll(".tickGroup")
			.data(y.domain(), d => d)
			.join(
				enter => {
					enter = enter.append("g")
						.attr("transform", yAxisTranslate)
						.attr("class", "tickGroup")
					enter.append("line")
						.attr("stroke", "#31363b")
						.attr("x2", 6)
					enter.append("g")
						.attr("class", "tickTextGroup")
						.call(pickLabel)
					return enter
				},
				update => update.call(
					update => {
						update.selectAll("g.tickTextGroup")
							.call(pickLabel)
						return update.transition(t)
							.attr("transform", yAxisTranslate)
					}),
				exit => exit.call(exit => exit.transition(100).remove())
			)
			.on("click", clickRemoveY)

		// Hide axes <path> elements
		heatmap.selectAll("path").style("opacity", "0");

		// Set up cell hit count toggle button.
		d3.select("#btn-toggle-counts")
			.on("click", () => {
				constants.showCounts = !constants.showCounts
				update(data)
			})

		// Hide/show genomic coordinates
		d3.select("#btn-toggle-scaffolds")
			.on("click", () => {
				constants.showScaffolds = !constants.showScaffolds
				update(data)
			})

		// Hide/show red borders for cells with multiple hits
		d3.select("#btn-toggle-count-borders")
			.on("click", () => {
				constants.showCountBorders = !constants.showCountBorders
				update(data)
			})

		// Hide/show red borders for cells with multiple hits
		d3.select("#btn-toggle-x-axis")
			.on("click", () => {
				constants.xAxisOnTop = !constants.xAxisOnTop
				update(data)
			})

		// Hide/show red borders for cells with multiple hits
		d3.select("#btn-toggle-multiline")
			.on("click", () => {
				constants.multiLineLabels = !constants.multiLineLabels
				update(data)
			})

		// Change width/height of heatmap cells
		d3.select("#input-cell-width")
			.on("change", function() {
				constants.cellWidth = +this.value;
				update(data)
			})
		d3.select("#input-cell-height")
			.on("change", function() {
				constants.cellHeight = +this.value;
				update(data)
			})

		// Wait until transition has finished until re-organising groups 
		setTimeout(() => {
			const offset = d3.max([40, heatmapX.node().getBBox().height]);
			heatmap.attr("transform", `translate(${constants.dendroWidth + 10}, ${offset})`)
			dendro.attr("transform", `translate(0, ${offset})`);
		}, 0)
	}

	update(data)
}
