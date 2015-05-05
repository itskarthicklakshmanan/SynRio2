$(document).ready(function() {
 
$(function() {
$( "#dialog" ).dialog({
autoOpen: false,
show: {
effect: "blind",
duration: 1000
},
hide: {
effect: "explode",
duration: 1000
}
});
$( "#mydiv" ).click(function() {
$( "#dialog" ).dialog( "open" );
});
});
 


 
  Shiny.addCustomMessageHandler("myCallbackHandler",
    function(color) {
      document.getElementById("res").innerHTML = color;
    }
  );

  document.getElementById("mydiv").onclick = function() {
    var number = Math.random();
    Shiny.onInputChange("mydata", number);
  };
 
 
});