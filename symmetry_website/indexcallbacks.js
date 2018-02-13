var submitData = function() {
    var synset = $(this).data('synset');
    var startIdx = $(this).data('start');
    var endIdx = $(this).data('end');
    $.get("writehtml.php",{category : synset, start: startIdx, end: endIdx}, function(data){
       console.log(data);
       if (data == "Done") {
           var table = document.getElementById("T"+synset).rows;
           var rowIdx = Math.floor(startIdx / 450);
           var row = table[rowIdx].cells;
           var cellIdx = Math.floor(startIdx / 45) - rowIdx*10;
           var cell = row[cellIdx];
           cell.children[0].style.background = "url('check_small.jpeg') no-repeat";
       } else {
           window.open(data, '_blank');
       }
    });
}

function checkButton(button) {
    var synset = $(button).data('synset');
    var startIdx = $(button).data('start');
    var endIdx = $(button).data('end');
    $.get("checkFile.php", {category : synset, start: startIdx, end: endIdx}, function(data){
        console.log(data);
        if (data == "exist") {
            console.log(button);
            button.style.background = "url('check_small.jpeg') no-repeat";
        }
    });
}

$(document).ready(function() {
    var buttons = document.getElementsByTagName("button");
    for (var i = 0; i < buttons.length; i++) {
        checkButton(buttons[i]);
    }
    $('button').click(submitData, submitData);
});

