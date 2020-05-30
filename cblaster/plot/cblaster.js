/* cblaster plot
 * Cameron L.M. Gilchrist, 2020.
 */

const constants = {
	"cellHeight": 27,
	"cellWidth": 44,
	"dendroWidth": 140,
}

if (typeof data === 'undefined') {
	// Load from HTTP server
	const data = d3.json("data.json").then(data => plot(data));
} else {
	plot(data);
}

// TODO: ignore zoom/drag transforms here
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
	const bbox = node.getBBox();
	node = node.cloneNode(true);
	node.setAttribute("width", bbox.width);
	node.setAttribute("height", bbox.height);
	node.setAttributeNS(xmlns, "xmlns", svgns);
	node.setAttributeNS(xmlns, "xmlns:xlink", xlinkns);
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
		.style("font-family", "sans-serif")
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
	return `<b>Query:</b> ${d.query}<br>`
		+ `<b>Organism:</b> ${data.labels[d.cluster].name}<br>`
		+ `<b>Scaffold:</b> ${data.labels[d.cluster].scaffold}<br>`
		+ `<b>${d.hits.length} ${d.hits.length > 1 ? "Hits:" : "Hit:"}</b>`
		+ '<table class="tooltip-table">'
		+ "<thead>"
		+ "<th>Name</th>"
		+ "<th>Ident. (%)</th>"
		+ "<th>Cov. (%)</th>"
		+ "<th>Bitscore</th>"
		+ "<th>E-value</th>"
		+ "</thead>"
		+ "<tbody>"
		+ d.hits.map(h => (
			"<tr>"
			+ `<td>${h.name}</td>`
			+ `<td>${h.identity.toFixed(2)}</td>`
			+ `<td>${h.coverage.toFixed(2)}</td>`
			+ `<td>${h.bitscore.toFixed(2)}</td>`
			+ `<td>${h.evalue}</td>`
			+ "</tr>"
		)).join("")
		+ "</tbody>"
		+ "</table>";
}

function tooltipOver(d) {
	// Heatmap cell mouseover behaviour to toggle tooltip opacity.
	d3.select(".tooltip").style("opacity", "1");
}

function tooltipMove(d, data) {
	// Heatmap cell mousemove behaviour to update content and screen position.
	d3.select(".tooltip")
		.html(getTooltipHTML(d, data))
		.style("left", d3.event.pageX + 10 + "px")
		.style("top", d3.event.pageY + 10 + "px")
}

function tooltipLeave(d) {
	// Heatmap cell mouseleave behaviour to toggle tooltip opacity.
	d3.select(".tooltip").style("opacity", "0");
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
	node.links().forEach(link => {
		link.source.y = y(link.source.data.length)	
		link.target.y = y(link.target.data.length || 0)	
	})
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

function addTooltipDiv(parent) {
	// Set up the <div> element used for the tooltip.
	const tooltip = parent.append("div")
		.classed("tooltip", true)
		.style("position", "absolute")
		.style("opacity", 0)
		.style("padding", "5px")
		.style("background-color", "white")
		.style("border", "solid")
		.style("border-width", "1px")
		.style("border-radius", "3px");
	return tooltip;
}

function pruneHierarchy(node, label) {
	/* Prune any children of a nested object whose name is equal to label.
	 * Additionally, if pruning results in an empty children list, remove the
	 * parent.
	*/
	return node.filter(child => {
		if (child.children) {
			child.children = pruneHierarchy(child.children, label)	
			if (child.children.length === 0)
				return false
		}
		return child.name !== label
	})
}

function plot(data) {
	data.matrix = flattenArray(data.matrix);

  const originalData = JSON.parse(JSON.stringify(data))

	// Reset to the original data. Have to make a deep copy here, since update
	// will mutate data
	d3.select("#btn-reset-filters")
		.on("click", () => {
			const copy = JSON.parse(JSON.stringify(originalData))
			update(copy)
		});

	const plotDiv = d3.select("#plot");
	const tooltip = addTooltipDiv(plotDiv);
	const svg = plotDiv.append("svg")
		.classed("wrapper-svg", true)
		.attr("id", "root_svg")
		.attr("cursor", "grab")
		.attr("xmlns", "http://www.w3.org/2000/svg");

	const g = svg.append("g").attr("transform", "translate(2,0)");

	svg.call(d3.zoom()
		.scaleExtent([0, 8])
		.on("zoom", () => g.attr("transform", d3.event.transform))
		.on("start", () => svg.attr("cursor", "grabbing"))
		.on("end", () => svg.attr("cursor", "grab"))
	);

	const x = d3.scaleBand().padding(0.05)
	const y = d3.scaleBand().padding(0.05)
	const colorScale = d3
		.scaleSequential(d3.interpolateBlues)
		.domain([0, 1])

	const def = svg.append("defs");
	def.append(() => colourBarGradient(colorScale));
	svg.append(() => colourBar());

	const dendro = g.append("g")
	const heatmap = g.append("g")
	const heatmapBG = heatmap.append("rect")
		.attr("class", "heatmap-bg")
		.style("fill", "white")
		.style("stroke", "black")
		.style("stroke-width", "thin");
	const cellGroup = heatmap.append("g")
	const heatmapX = heatmap.append("g")
		.attr("class", "g-heatmap-xaxis")
	const heatmapY = heatmap.append("g")
		.attr("class", "g-heatmap-yaxis")

	// Save SVG
	d3.select("#btn-save-svg")
		.on("click", () => {
			const blob = serialise(svg)
			download(blob, "cblaster.svg")
		})

	const summaryHTML = getSummaryHTML(data.counts)
	d3.select("#p-result-summary")
		.html(summaryHTML)

	let showCounts = true;
	d3.select("#btn-toggle-counts")
		.on("click", () => {
			showCounts = !showCounts
			let result = showCounts ? "visible" : "hidden"
			d3.selectAll(".heatmap-count")
				.style("visibility", result)
		})

	function update(data) {
		t = d3.transition().duration(500)

		x.domain(data.queries)
			.range([0, data.queries.length * constants.cellWidth])

		const tree = getScaledTree(
			data.hierarchy,
			constants.dendroWidth,
			constants.cellHeight * Object.keys(data.labels).length,
		);

		dendro.selectAll("path")
			.data(tree.links())
			.join("path")
			.attr("fill", "none")
			.attr("stroke", "#555")
			.attr("stroke-opacity", 0.4)
			.attr("stroke-width", 1.5)
			.attr("d", elbow);

		const dendroBBox = dendro.node().getBBox()
		const dendroLeaves = tree.leaves()

		y.domain(dendroLeaves.map(l => l.data.name))
			.range([0, dendroLeaves.length * constants.cellHeight])

		const translate = (d) => `translate(${x(d.query)}, ${y(d.cluster)})`

		const cells = cellGroup
			.selectAll("g")
			.data(data.matrix, d => `${d.query}-${d.cluster}`)
			.join(
				enter => {
					enter = enter.append("g").attr("transform", translate)
					enter.append("rect")
					enter.append("text")
					return enter
				},
				update => update.call(
					update => update.transition(t)
						.attr("transform", translate)
				)
			)
			.filter(d => d.hits.length > 0)
			.on("mouseover", tooltipOver)
			.on("mousemove", d => tooltipMove(d, data))
			.on("mouseleave", tooltipLeave)

		cells.selectAll("rect")
			.attr("fill", d => colorScale(d.value / 100))
			.attr("width", x.bandwidth())
			.attr("height", y.bandwidth())
			.style("stroke-width", "thin")
			.style("stroke", d => d.hits.length > 1 ? "red" : null)

		cells.selectAll("text")
			.text(d => d.hits.length)
			.attr("x", constants.cellWidth / 2)
			.attr("y", constants.cellHeight / 2 + 4)
			.attr("class", "heatmap-count")
			.style("text-anchor", "middle")
			.style("fill", d => d.value < 60 ? colorScale(1) : colorScale(0))

		const bbox = cellGroup.node().getBBox()
		const heatWidth = data.queries.length * constants.cellWidth
		const heatHeight = Object.keys(data.labels).length * constants.cellHeight

		heatmapBG.transition(t)
			.attr("width", heatWidth)
			.attr("height", heatHeight)

		const axisTop = d3.axisTop(x)
			.tickSize(6)
			.tickPadding(6)

		heatmapX.transition(t)
			.call(axisTop)

		heatmapX.selectAll("text")
			.style("font-size", "12px")
			.style("text-anchor", "start")
			.attr("transform", "rotate(-25)")
			.on("click", (query) => {
				let newData = {
					...data,
					"queries": data.queries.filter(q => q !== query),
					"matrix": data.matrix.filter(cell => cell.query !== query)
				}
				update(newData)
			})

		const axisRight = d3.axisRight(y)
			.tickSize(6)

		heatmapY.transition(t)
			.attr("transform", `translate(${heatWidth}, 0)`)
	
		heatmapY.call(axisRight)
			.selectAll("text")
				.call(t => {
					t.text("")
					t.append("tspan")
						.text(d => data.labels[d].name)
						.attr("x", 10)
						.attr("dy", "-0.1em")
						.attr("text-anchor", "start")
						.style("font-family", "sans-serif")
						.style("font-style", "italic")
						.style("font-size", "12px")
					t.append("tspan")
						.text(d => data.labels[d].scaffold)
						.attr("x", 10)
						.attr("dy", "1.3em")
						.attr("text-anchor", "start")
						.style("font-style", "normal")
						.style("font-family", "sans-serif")
						.style("font-size", "10px")
						.style("fill", "grey");				
				})
				.on("click", (label) => {
					const newData = {
						...data,
						"hierarchy": {
							...data.hierarchy,
							"children": pruneHierarchy(data.hierarchy.children, label)
						},
						"matrix": data.matrix.filter(cell => cell.cluster !== label)
					}
					delete newData.labels[label]
					update(newData)
				})

		heatmap.selectAll("path").style("opacity", "0");

		// Wait until transition has finished until re-organising groups 
		setTimeout(() => {
			const offset = d3.max([40, heatmapX.node().getBBox().height]);
			heatmap.attr("transform", `translate(${constants.dendroWidth + 10}, ${offset})`)
			dendro.attr("transform", `translate(0, ${offset})`);
		}, 0)
	}

	update(data)
}
