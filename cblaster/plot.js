/* cblaster plot
 *
 * Expects JSON data in form:
 *
 * dendrogram = [
 *   {
 *     "name": "cluster1_key",
 *     "children": [
 *       { "name": cluster2_key", "length": 0.8 },
 *       { "name": cluster3_key", "length": 1.2 },
 *     ]
 *   },
 *   { ... },
 * ]
 *
 * clusters = [
 *   {
 *     "key": "cluster1_key",
 *     "query_1": [
 *       { "name": "hit1_name", "identity": 89.23 },
 *       { "name": "hit2_name", "identity": 72.58 },
 *     ],
 *     "query_2": [ ... ],
 *   },
 *   { "key": "cluster2_key", ... },
 * ]
 *
 * Cameron L.M. Gilchrist, 2020.
 */

const data = d3
	.json("plot.json")
	.then(data => {
		plot(data);
	});


/* Scale the branch lengths of the cluster dendrogram.
 * This requires cluster distances to be normalised (e.g. divide all distances
 * in SciPy linkage matrix by the maximum distance). A linear scale maps the
 * normalised distances to the total width of the dendrogram, and then y values
 * of every link source and target are scaled accordingly. Terminal targets will
 * have no length attribute and thus default to 0 (i.e. drawn to max width).
 */
function scaleBranchLengths(node, w) {
	const y = d3.scaleLinear()
		.domain([1, 0])  // inversed, since thats how Scipy reports distances
		.range([0, w]);
	node.links().map(link => {
		link.source.y = y(link.source.data.length)	
		link.target.y = y(link.target.data.length || 0)	
	})
}

function plot(data) {

	// Make a copy of the data to restore from when elements are hidden.
	const originalData = data;

	const cellHeight = 27;
	const cellWidth = 44;
	const width = 1000;
	const height = 1000;

	const queries = data.queries;
	const clusters = data.matrix.map(d => d.name);
	const scaffolds = data.matrix.map(d => d.scaffold);

	const dendroWidth = 120;
	const dendroHeight = clusters.length * cellHeight;

	// x-axis scale is only used for the heatmap, and is based on user queries.
	const x = d3.scaleBand()
		.domain(queries)
		.range([0, queries.length * cellWidth])
		.padding(0.05);

	// y-axis scale is used for both plots, and is based on hit clusters.
	const y = d3.scaleBand()
		.domain(clusters)
		.range([0, dendroHeight])
		.padding(0.05);

	// Heat map color scheme
	const colors = d3
		.scaleSequential(d3.interpolateBlues)
		.domain([0, 1]);

	const elbow = (d) => {
		return `
			M${d.source.y},${d.source.x}
			V${d.target.x}
			H${d.target.y}
		`;
	}

	// Generate the dendrogram from the clustering hierarchy.
	const tree = (data) => {
		const root = d3.hierarchy(data);
		return d3.cluster()
			.size([dendroHeight, dendroWidth])
			.separation(() => 1)(root);
	}

	const root = tree(data.hierarchy);
	scaleBranchLengths(root, dendroWidth);

	const svg = d3
		.select("#plot")
		.append("svg")
		.attr("id", "root_svg")
		.attr("width", width)
		.attr("height", height);

	// Define a linear gradient for the color bar
	const defs = svg.append("defs");
	const cbarGradient = defs.append("linearGradient")
		.attr("id", "cbar-gradient")
		.attr("x1", "0%")
		.attr("x2", "100%");
	cbarGradient.append("stop")
		.attr("offset", "0%")
		.attr("stop-color", "white");
	cbarGradient.append("stop")
		.attr("offset", "100%")
		.attr("stop-color", colors(1));

	// Create the color bar, as well as text labels
	const cbar = svg.append("g")
		.attr("transform", "translate(2, 2)");
	cbar
		.append("rect")
		.attr("width", 150)
		.attr("height", 12)
		.style("fill", "url(#cbar-gradient)");
	cbar
		.append("rect")
		.style("fill", "none")
		.style("outline", "thin solid black")
		.attr("width", 150)
		.attr("height", 10)
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
	cbar
		.selectAll("text")
		.style("font-family", "sans-serif")
		.style("font-size", "12px")

	const g = svg.append("g").attr("transform", "translate(2,0)");

	// svg.call(d3.zoom()
	// 	.extent([[0, 0], [width, height]])
	// 	.scaleExtent([1, 8])
	// 	.on("zoom", () => g.attr("transform", d3.event.transform))
	// )

	// Draw the paths in the dendrogram
	const dendro = g.append("g")
		.attr("transform", "translate(0, 50)")
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

	// Text labels for the dendrogram
	dendro.append("g")
    .selectAll("text")
    .data(root.leaves())
		.join("g")
	.append("text")
		.attr("x", d => d.y)
		.attr("y", d => d.x)
		.attr("dy", "0em")
		.attr("dx", ".5em")
		.style("font-family", "sans-serif")
		.style("font-style", "italic")
		.style("font-size", "12px")
    .filter(d => !d.children)
		.text(d => d.data.name)
		.attr("text-anchor", "start")
	.append("svg:tspan")
		.text(d => d.data.scaffold)
		.style("font-style", "normal")
		.style("font-family", "sans-serif")
		.style("font-size", "10px")
		.style("fill", "grey")
		.attr("x", d => d.y)
		.attr("y", d => d.x)
		.attr("dy", "1.2em")
		.attr("dx", ".5em")
		.attr("text-anchor", "start");

	// Draw the heatmap
	const cells = [];
	data.matrix.forEach((row, i) => {
		row.values.forEach((query, i) => {
			cells.push({
				query: query.name,
				cluster: row.name,
				hits: query.hits,
				identity: Math.max(...query.hits.map(h => h.identity))
			});
		})
	});

	const tooltip = d3
		.select("#plot")
		.append("div")
		.attr("class", "tooltip")
		.style("position", "absolute")
		.style("opacity", 0)
		.style("padding", "5px")
		.style("background-color", "white")
		.style("border", "solid")
		.style("border-width", "1px")
		.style("border-radius", "3px");

	const mouseover = (d) => {
		tooltip.style("opacity", "1");
	};

	const mousemove = (d) => {
		tooltip
			.html(
				`<b>${d.hits.length} ${d.hits.length > 1 ? "Hits:" : "Hit:"}</b>`
				+ '<ol style="list-style: none; margin: 0; padding-left: .3em;">'
				+ d.hits.map(h => '<li style="font-family: sans-serif; font-size: 12px;">' + h.name + ", " + (h.identity * 100).toFixed(2) + "% identity</li>").join("")
				+ "</ol>"
			)
			.style("left", d3.event.pageX + 10 + "px")
			.style("top", d3.event.pageY + 10 + "px")
	};

	const mouseleave = (d) => {
		tooltip.style("opacity", "0");
	};

	const dendroBBox = d3.select(".g-dendro").node().getBBox()
	const heatmap = g.append("g")
		.attr("transform", `translate(${dendroBBox.x + dendroBBox.width + 10}, 50)`);

  const heatGroup = heatmap.append("g")
		.classed("g-heatmap", true);

	heatGroup
		.selectAll("rect")
		.data(cells)
		.enter()
		.append("rect")
		.attr("fill", d => colors(d.identity))
		.attr("width", x.bandwidth())
		.attr("height", y.bandwidth())
		.attr("x", d => x(d.query))
		.attr("y", d => y(d.cluster))
		.on("mouseover", mouseover)
		.on("mousemove", mousemove)
		.on("mouseleave", mouseleave);

	const heatBBox = d3.select(".g-heatmap").node().getBBox()
	heatmap.append("rect")
		.attr("width", heatBBox.width + 4)
		.attr("height", heatBBox.height + 4)
		.style("fill", "none")
		.style("outline", "thin solid black");

	const x_axis = d3
		.axisTop(x)
		.tickSize(0)
		.tickPadding(6);

	const y_axis = d3
		.axisLeft(y)
		.tickSize(0)
		.tickFormat((d, i) => scaffolds[i])
		.tickPadding(6);

	heatmap
		.append("g")
		.call(x_axis)
		.selectAll("text")
		.style("font-size", "12px")
		.style("text-anchor", "start")
		.attr("transform", "rotate(-25)");

	// Hide domain <path> element created in axisTop()
	heatmap.selectAll("path")
		.style("opacity", "0");
}
