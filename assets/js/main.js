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

