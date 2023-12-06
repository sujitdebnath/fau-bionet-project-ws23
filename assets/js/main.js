document.addEventListener("DOMContentLoaded", function(event) {
	// document.getElementById('run').onclick = function() {
    //     	alert("button was clicked");
    //     	document.getElementById("run").disabled = true;
	// }
});

function freeze_buttons() {
	// Disables the 'run' button so user cannot push it.
	document.getElementById("run").disabled = true;
	document.getElementById("plot_figures").disabled = true;
}

function release_buttons() {
	// Disables the 'run' button so user cannot push it.
	document.getElementById("run").disabled = false;
	document.getElementById("plot_figures").disabled = false;
}

function auto_plot() {
	document.getElementById("plot_figures").click();
}

function goto_first_tab() {
	var navLinks = document.querySelectorAll('a.nav-link');

	// Now you can loop through the navLinks NodeList and perform actions on each link
	navLinks.forEach(function(link) {
		// For example, you can log the text content of each link
		if (link.textContent === "Highest Expressed Genes") {
			link.click();
		}
	});
}