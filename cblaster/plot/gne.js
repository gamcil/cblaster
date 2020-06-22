/* Gene neighbourhood plot.
 * */

if (typeof data === 'undefined') {
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
	const bbox = node.getBBox();
	node = node.cloneNode(true);
	d3.select(node).select(".marker").style("opacity", 0)
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

function plot(data) {
	const width = 900
	const height = 600
	const margin = {
		left: 160,
		right: 30,
		bottom: 60,
		top: 20,
	}

	// Add <svg> element. Viewport has height of 100vh such that the SVG will
	// always expand, at the viewBox aspect ratio, to fill the screen vertically
	const svg = d3.select("#plot")
		.append("svg")
		.attr("xmlns", "http://www.w3.org/2000/svg")
		.attr("height", "100vh")
		.attr("viewBox", [0, 0, width, height])

	// Set up the save SVG button
	d3.select("#btn-save-svg")
		.on("click", () => {
			const blob = serialise(svg)
			download(blob, "gne.svg")
		})

	// Add a <rect> element to use as a clip path for the chart, since zoom/pan
	// transformations will draw the lines outside of the chart area.
	const clip = svg.append("clipPath")
		.attr("id", "clip")
		.append("rect")
		.attr("x", margin.left)
		.attr("y", margin.top)
		.attr("width", width - margin.left - margin.right)
		.attr("height", height - margin.top - margin.bottom)

	// Chart skeleton
	const g = svg.append("g")
	const chart = g.append("g")		
	const yAxes = chart.append("g")
	const clusters = chart.append("g")
	const means = chart.append("g")
	const medians = chart.append("g")

	// Need min/max of data points for scale domains
	const [minGap, maxGap] = d3.extent(data, d => d.gap)
	const [minMean, maxMean] = d3.extent(data, d => d.means)
	const [minClusters, maxClusters] = d3.extent(data, d => d.clusters)

	// TODO: add linear/log scale toggle
	//    If doing this, start from first non-zero mean/median value (as below)
	//		Then, make sure not possible for d.gap to be 0 when using the log scale
	//		Will have to adjust x values of lines, axis, whatever is dependent
	// let blah = 0
	// for (let i = 0; i < data.length; i++) {
	// 	if (data[i].gap > blah) {
	// 		blah = data[i].gap
	// 		break
	// 	}
	// }

	// x scale from intergenic distance samples
	const x = d3.scaleLinear()
		.domain([minGap, maxGap])
		.range([margin.left, width - margin.right])
	const xAxis = (g, x) => g
		.attr("transform", `translate(0, ${height - margin.bottom})`)
		.call(d3.axisBottom(x))
	const gx = chart.append("g")
		.call(xAxis, x)

	const top = (height) => height * 0.7 - margin.bottom

	// y scale for mean/median cluster size
	const y0 = d3.scaleLinear()
		.domain([minMean, maxMean])
		.range([top(height) - 20, margin.top])
		.nice()
	const yAxis0 = yAxes.append("g")
		.attr("transform", `translate(${margin.left}, 0)`)
		.call(d3.axisLeft(y0))

	// y scale for no. clusters
	const y1 = d3.scaleLinear()
		.domain([minClusters, maxClusters])
		.range([height - margin.bottom - 20, top(height)])
		.nice()
	const yAxis1 = yAxes.append("g")
		.attr("transform", `translate(${margin.left}, 0)`)
		.call(d3.axisLeft(y1))

	// Add chart-width horizontal bars
	yAxis0.append("g")
		.selectAll("line")
		.data(y0.ticks())
		.join("line")
		.attr("class", "gridline")
		.attr("stroke", "#ddd")
		.attr("x1", 0)
		.attr("x2", width - margin.left - margin.right)
		.attr("y1", d => y0(d))
		.attr("y2", d => y0(d))
	yAxis1.append("g")
		.selectAll("line")
		.data(y1.ticks())
		.join("line")
		.attr("class", "gridline")
		.attr("stroke", "#ddd")
		.attr("x1", 0)
		.attr("x2", width - margin.left - margin.right)
		.attr("y1", d => y1(d))
		.attr("y2", d => y1(d))
	yAxis0.selectAll("line:not(.gridline), path, text")
		.attr("transform", "translate(-6, 0)")
		.raise()
	yAxis1.selectAll("line:not(.gridline), path, text")
		.attr("transform", "translate(-6, 0)")
		.raise()

	// <line> d attribute generators
	const cLineD = (data, x) => d3.line()
		.x(d => x(d.gap))
		.y(d => y1(d.clusters))
		(data)
	const meanLineD = (data, x) => d3.line()
		.x(d => x(d.gap))
		.y(d => y0(d.means))
		(data)
	const mediansLineD = (data, x) => d3.line()
		.x(d => x(d.gap))
		.y(d => y0(d.medians))
		(data)

	// Draw lines for mean/median cluster size and total clusters
	const clusterLine = clusters.append("path")
		.datum(data)
		.attr("clip-path", "url(#clip)")
		.attr("fill", "none")
		.attr("stroke", "red")
		.attr("d", cLineD(data, x))
	const meanLine = means.append("path")
		.datum(data)
		.attr("clip-path", "url(#clip)")
		.attr("fill", "none")
		.attr("stroke", "green")
		.attr("d", meanLineD(data, x))
	const medianLine = medians.append("path")
		.datum(data)
		.attr("clip-path", "url(#clip)")
		.attr("fill", "none")
		.attr("stroke", "blue")
		.attr("d", mediansLineD(data, x))

	// Set up hovering tooltip/marker line. Uses d3.bisector to find nearest data
	// index to a given mouse position based on intergenic gap value.
	const bisect = d3.bisector(d => d.gap)
	const tooltip = d3.select("#plot")
		.append("div")
		.style("background", "rgba(0, 0, 0, 0.1)")
		.style("opacity", 0)
		.style("padding", "5px")
		.style("position", "absolute")
	const marker = chart.append("line")
		.style("stroke", "#999")
		.style("stroke-width", "0.5")
		.style("stroke-dasharray", "5 3")
		.attr("class", "marker")
		.attr("y1", margin.top)
		.attr("y2", height - margin.bottom)
		.attr("x1", -100)
		.attr("x2", -100)

	// Transformed x axis; altered on zoom/pan
	let xz = x

	// Attach tooltip behaviour
	svg.on("mousemove click touchmove", function() {
		const m = d3.mouse(this)
		if (m[0] > width - margin.right) return
		const gap = xz.invert(m[0])
		const i = bisect.right(data, gap)
		if (i < data.length) {
			marker
				.attr("x1", xz(data[i].gap))
				.attr("x2", xz(data[i].gap))
			let bbox = marker.node().getBoundingClientRect()
			tooltip
				.html(`
				<table>
					<tbody>
						<tr>
							<td><b>Gap size</b></td>
							<td>${data[i].gap}</td>
						</tr>
						<tr style="color: green">
							<td><b>Mean</b></td>
							<td>${parseInt(data[i].means)}</td>
						</tr>
						<tr style="color: blue">
							<td><b>Median</b></td>
							<td>${data[i].medians}</td>
						</tr>
						<tr style="color: red">
							<td><b>Clusters</b></td>
							<td>${data[i].clusters}</td>
						</tr>
					</tbody>
				</table>
				`)
			tooltip
				.style("opacity", 1)
				.style("left", bbox.left + 5 + "px")
				.style("top", bbox.height / 2 - bbox.y + "px")
		}
	})

	// Define and attach zoom/pan behaviour. Rescales x axis, then updates each
	// line and the x axis group.
	const zoomed = () => {
		xz = d3.event.transform.rescaleX(x)
		gx.call(xAxis, xz)
		meanLine.attr("d", meanLineD(data, xz))
		medianLine.attr("d", mediansLineD(data, xz))
		clusterLine.attr("d", cLineD(data, xz))
	}
	const zoom = d3.zoom()
		.scaleExtent([1, 32])
		.translateExtent([[margin.left, 0], [width - margin.right, height]])
		.extent([[margin.left, -Infinity], [width - margin.right, Infinity]])
		.on("zoom", zoomed)
	svg.call(zoom)

	// Add chart axis titles
	chart.append("text")
		.text("Intergenic distance threshold (bp)")
		.attr("x", (margin.left + width) / 2)
		.attr("y", margin.top + height - margin.bottom / 2)
		.attr("text-anchor", "middle")
		.style("font-size", "12px")
	chart.append("text")
		.text("Cluster size (bp)")
		.attr("x", margin.left - 50)
		.attr("y", () => {
			let [a, b] = y0.range()
			return (a - b) / 2
		})
		.attr("fill", "black")
		.attr("text-anchor", "end")
		.style("font-size", "12px")
	chart.append("text")
		.text("Total clusters")
		.attr("x", margin.left - 50)
		.attr("y", () => {
			let [a, b] = y1.range()
			return top(height) + (a - b) / 2
		})
		.attr("fill", "black")
		.attr("text-anchor", "end")
		.style("font-size", "12px")
}
