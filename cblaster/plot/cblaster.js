/* cblaster plot
 * Cameron L.M. Gilchrist, 2020.
 */


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
	console.log(svg);
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

function resetFilters(originalData) {
	/* Resets any hidden rows/columns to original data.
	 * This is attached to the Reset filters button in plot.html.
	*/
	d3.select("#plot")
		.select("svg")
		.remove();
	plot(originalData)
}

function scaleBranchLengths(node, w) {
	/* Scales the branch lengths of the cluster dendrogram.
	 * This requires cluster distances to be normalised (e.g. divide all distances
	 * in SciPy linkage matrix by the maximum distance). A linear scale maps the
	 * normalised distances to the total width of the dendrogram, and then y values
	 * of every link source and target are scaled accordingly. Terminal targets will
	 * have no length attribute and thus default to 0 (i.e. drawn to max width).
	 */
	const y = d3.scaleLinear()
		.domain([1, 0])  // inversed, since thats how Scipy reports distances
		.range([0, w]);
	node.links().map(link => {
		link.source.y = y(link.source.data.length)	
		link.target.y = y(link.target.data.length || 0)	
	})
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

function plot(data) {
	// Save data being plotted, and bind to button to allow resetting of any
	// hidden rows/columns, transforms, etc
	const originalData = data;
	d3.select("#btn-reset-filters")
		.on("click", () => resetFilters(originalData));

	// Define constant size for heatmap cells. Scales are based on these values
	// * the number of total clusters/queries, so the figure will not scale to fit
	// a given window size.
	const cellHeight = 27;
	const cellWidth = 44;

	const queries = data.queries;
	const clusters = data.matrix.map((d, i) => i);

	// Define width and height of the dendrogram.
	const dendroWidth = 140;
	const dendroHeight = clusters.length * cellHeight;

	// x-axis scale is only used for the heatmap, and is based on user queries.
	const x = d3.scaleBand()
		.domain(data.queries)
		.range([0, queries.length * cellWidth])
		.padding(0.05);

	// Heat map color scheme
	const colors = d3
		.scaleSequential(d3.interpolateBlues)
		.domain([0, 1]);

	// Create the root SVG element
	const svg = d3
		.select("#plot")
		.append("svg")
		.classed("wrapper-svg", true)
		.attr("id", "root_svg")
		.attr("cursor", "grab")
		.attr("xmlns", "http://www.w3.org/2000/svg");

	// Create a group which everything in the plot will be inside
	const g = svg.append("g").attr("transform", "translate(2,0)");

	// Set up pan/zoom behaviour. Affects the <g> element, which does not include
	// the colour bar.
	svg.call(d3.zoom()
		.scaleExtent([0, 8])
		.on("zoom", () => g.attr("transform", d3.event.transform))
		.on("start", () => svg.attr("cursor", "grabbing"))
		.on("end", () => svg.attr("cursor", "grab"))
	);

	// Create the colour bar
	const def = svg.append("defs");
	def.append(() => colourBarGradient(colors));
	svg.append(() => colourBar());

	// Generate the dendrogram from the clustering hierarchy, then scale branch
	// lengths based on distances from the linkage matrix.
	const tree = (data) => {
		const root = d3.hierarchy(data);
		return d3.cluster()
			.size([dendroHeight, dendroWidth])
			.separation(() => 1)(root);
	}
	const root = tree(data.hierarchy);
	scaleBranchLengths(root, dendroWidth);

	// y-axis scale is used for both plots, and is based on hit clusters.
	// Must be calculated after creating the tree, so we can get order.
	const y = d3.scaleBand()
		.domain(root.leaves().map(l => l.data.name))
		.range([0, dendroHeight])
		.padding(0.05);

	// Draw the paths in the dendrogram
	const dendro = g.append("g")
		.classed("g-dendro", true);
	dendro.append("g")
    .attr("fill", "none")
    .attr("stroke", "#555")
    .attr("stroke-opacity", 0.4)
    .attr("stroke-width", 1.5)
		.selectAll("path")
    .data(root.links())
    .join("path")
		.attr("d", elbow);

	// Generate the heatmap cell tooltip, hidden by default.
	const tooltip = d3
		.select("#plot")
		.append("div")
		.classed("tooltip", true)
		.style("position", "absolute")
		.style("opacity", 0)
		.style("padding", "5px")
		.style("background-color", "white")
		.style("border", "solid")
		.style("border-width", "1px")
		.style("border-radius", "3px");

	const dendroBBox = d3.select(".g-dendro").node().getBBox()

	// Define groups for the heatmap.
	// One for entire heatmap, another just for heatmap cells.
	const heatmap = g.append("g");

	const heatWidth = cellWidth * data.queries.length;
	const heatHeight = cellHeight * data.labels.length + 2;

	// Draw a padded border around the heatmap
	heatmap.append("rect")
		.attr("width", heatWidth)
		.attr("height", heatHeight)
		.style("fill", "white")
		.style("stroke", "black")
		.style("stroke-width", "thin");

  const heatGroup = heatmap.append("g");
	const heatCell = heatGroup.selectAll("rect")
		.data(flattenArray(data.matrix))
		.enter()
		.append("g");

	// Generates the HTML content for a cell hovering tooltip.
	// It provides the name of the query, as well as a table of each hit.
	const getTooltipHTML = (d) => {
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

	// Heatmap cell mouseover behaviour to toggle tooltip opacity.
	const tooltipOver = (d) => {
		tooltip.style("opacity", "1");
	}

	// Heatmap cell mousemove behaviour to update content and screen position.
	const tooltipMove = (d) => {
		const content = getTooltipHTML(d);
		tooltip.html(content)
			.style("left", d3.event.pageX + 10 + "px")
			.style("top", d3.event.pageY + 10 + "px")
	}

	// Heatmap cell mouseleave behaviour to toggle tooltip opacity.
	const tooltipLeave = (d) => {
		d3.select(".tooltip").style("opacity", "0");
	}

	heatCell.filter(d => d.hits.length > 0)
		.on("mouseover", tooltipOver)
		.on("mousemove", tooltipMove)
		.on("mouseleave", tooltipLeave);

	// Draw heatmap cells
	heatCell.filter(d => d.hits.length > 0)
		.append("rect")
		.classed("heatmap-cell", true)
		.attr("fill", d => colors(d.value / 100))
		.attr("width", x.bandwidth())
		.attr("height", y.bandwidth())
		.attr("x", d => x(d.query))
		.attr("y", d => y(d.cluster))
		.style("stroke-width", "thin")
		.style("stroke", d => d.hits.length > 1 ? "red" : null);

	// Draw hit counts in the center of each cell
	heatCell.filter(d => d.hits.length > 0)
		.append("text")
		.classed("heatmap-count", true)
		.text(d => d.hits.length)
		.attr("x", d => x(d.query) + cellWidth / 2)
		.attr("y", d => y(d.cluster) + cellHeight / 2 + 4)
		.style("text-anchor", "middle")
		.style("fill", d => d.value < 60 ? colors(1) : colors(0))

	const heatBBox = heatGroup.node().getBBox()

	// Species names and scaffold coordinates
	heatmap.append("g")
    .selectAll("text")
    .data(root.leaves())
		.join("g")
	.append("text")
		.text(d => data.labels[d.data.name].name)
		.attr("x", heatBBox.x + heatWidth)
		.attr("y", d => d.x)
		.attr("dy", "0em")
		.attr("dx", "1.1em")
		.attr("text-anchor", "start")
		.style("font-family", "sans-serif")
		.style("font-style", "italic")
		.style("font-size", "12px")
	.append("svg:tspan")
		.text(d => data.labels[d.data.name].scaffold)
		.attr("x", heatBBox.x + heatWidth)
		.attr("y", d => d.x)
		.attr("dy", "1.2em")
		.attr("dx", "1.2em")
		.attr("text-anchor", "start")
		.style("font-style", "normal")
		.style("font-family", "sans-serif")
		.style("font-size", "10px")
		.style("fill", "grey");

	// Set up ticks and labels on top of the heatmap. The group is defined
	// separately to the axisTop call since we use its bounding box to determine
	// the offset from the top of the <svg>.
	const heatXAxis = heatmap.append("g");
	heatXAxis.call(d3.axisTop(x).tickSize(6).tickPadding(6))
		.selectAll("text")
		.style("font-size", "12px")
		.style("text-anchor", "start")
		.attr("transform", "rotate(-25)");

	// Set up ticks on the right of the heatmap, and remove any labels since
	// those are drawn separately.
	const heatYAxis = heatmap.append("g")
		.call(d3.axisRight(y).tickSize(6))
		.attr("transform", `translate(${heatWidth}, 0)`)
		.selectAll("text")
		.remove();

	// Get height of heatmap ticks, then translate dendrogram/heatmap y values accordingly.
	const offset = d3.max([40, heatXAxis.node().getBBox().height]);
	heatmap.attr("transform", `translate(${dendroBBox.x + dendroBBox.width + 6}, ${offset})`)
	dendro.attr("transform", `translate(0, ${offset})`);

	// Hide domain <path> element created in axisTop()
	heatmap.selectAll("path").style("opacity", "0");

	// Bind save button to serialisation/download
	d3.select("#btn-save-svg")
		.on("click", () => {
			const blob = serialise(svg.node());
			download(blob, "cblaster.svg")
		})

	// Generate HTML for the count summary
	d3.select("#p-result-summary").html(() => getSummaryHTML(data.counts));

	// Toggle visibility of cell counts
	let showCounts = true;
	d3.select("#btn-toggle-counts")
		.on("click", () => {
			showCounts = !showCounts;
			let result = showCounts ? "visible" : "hidden";
			d3.selectAll(".heatmap-count")
				.style("visibility", result);
			})
}
