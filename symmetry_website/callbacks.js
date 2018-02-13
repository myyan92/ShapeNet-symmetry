var states = new Array(45);    // 0:correct, 1:question, 2:wrong
for (var i = 0; i < 45; i++) {states[i]=0;}
var sourceSwap = function() {
    var newSource = $(this).data('alt-src');
    $(this).data('alt-src', $(this).attr('src'));
    $(this).attr('src', newSource);
}

var click = function() {
    var cell = this.parentNode.parentNode.parentNode.parentNode.parentNode.parentNode;
    var column = cell.cellIndex;
    var row = cell.parentNode.rowIndex;
    var icon = cell.children[0].rows[0].cells[0].children[0].children[1];
    var currentName = icon.src;
    console.log(currentName);
    if (currentName.endsWith("question.jpeg")) {
        icon.src = "cross.jpeg";
        states[row*3+column] = 2;
    } else if (currentName.endsWith("cross.jpeg")) {
        icon.src = "check.jpeg";
        states[row*3+column] = 0;
    } else {
        icon.src = "question.jpeg";
        states[row*3+column] = 1;
    }
}

var submitData = function() {
    var synset = $(this).data('synset');
    var startIdx = $(this).data('start');
    var endIdx = $(this).data('end');
    console.log(synset);
    console.log(startIdx);
    console.log(endIdx);
    $.get("logData.php",{category : synset, start: startIdx, end: endIdx, annotation : states}, function(data){
       console.log(data);
       window.open(data, '_self');
    });
}

window.onload = function() {
    console.log("loading hover images");
    imgs = document.getElementsByTagName("img");
    for (var i = 0; i < imgs.length; i++) {
        preloadImage = new Image();
        preloadImage.src = imgs[i].getAttribute("data-alt-src");
    }
}

$(document).ready(function() {
    $('img.withHover').hover(sourceSwap, sourceSwap);
    $('img').click(click, click);
    $('button').click(submitData, submitData);
});
