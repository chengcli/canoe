$(document).ready(function() {
  $("ol.publist li").each(function() {
    $(this).html(function(index, text) {
      // Cheng Li
      text = text.replace('Cheng Li', '<a href="https://chengcli.io">Cheng Li</a>');
      text = text.replace('C. Li', '<a href="https://chengcli.io">C. Li</a>');
      // Huazhi Ge
      text = text.replace('Huazhi Ge', '<a href="http://people.ucsc.edu/~hge2">Huazhi Ge</a>');
      // Xi Zhang
      text = text.replace('X. Zhang', '<a href="https://ucscplanetaryscience.com/xi-zhang/">X. Zhang</a>');
      // Tianhao Le
      text = text.replace('T. Le', '<a href="https://happysky19.github.io/">T. Le</a>');
      return text;
    });
  });
});
