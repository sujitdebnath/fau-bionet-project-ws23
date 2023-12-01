document.addEventListener("DOMContentLoaded", function(event) {
	document.getElementById('run').onclick = function() {
        	alert("button was clicked");
        	document.getElementById("run").disabled = true;
	}	
});

