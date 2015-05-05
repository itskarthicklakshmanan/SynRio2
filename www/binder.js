
$(document).ready(function() {

setInterval(function(){
  if ($('html').attr('class')=='shiny-busy') {
    setTimeout(function() {
      if ($('html').attr('class')=='shiny-busy') {
        $('div.busy').show()
      }
    }, 100)
  } else {
    $('div.busy').hide()
  }
}, 100)

  setTimeout(function () {
      $('div.intro').hide();
  }, 17000);
  
    setTimeout(function () {
      $('div.intro1').hide();
  }, 8000);

 
});