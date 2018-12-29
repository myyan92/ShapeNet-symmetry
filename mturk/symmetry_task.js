// Helper functions
function swap_image_source() {
    var newSource = $(this).data('alt-src');
    $(this).data('alt-src', $(this).attr('src'));
    $(this).attr('src', newSource);
}

function preloadImages(images) {
  // Preload all images
  var imgs = [];
  _.each(images, function(img_url) {
    var img = new Image();
    img.onload = function() { console.log('loaded image from ' + img_url); };
    img.src = img_url;
  });
}

function preloadAllImages(items) {
  preloadImages(items.map(function(x) { return x.image1; }));
  preloadImages(items.map(function(x) { return x.image2; }));  
}

function populateCategoryTable(categoriesTable, categories, items) {
  var itemsByExpected = _.groupBy(items, function(x) { return x.expected; });
  for (var i = 0; i < categories.length; i++) {
    var category = categories[i];
    var row = $('<tr></tr>');
    row.append($('<td></td>').text(category.label));
    var details = $('<td></td');
    details.append($('<span><span>').text(category.description));
    row.append(details);
    var rowItems = itemsByExpected[category.value];
    // console.log('items for ', category, rowItems);
    if (rowItems) {
      details.append('<br/>');
      for (var j = 0; j < rowItems.length; j++) {
        var img = createItemImage(rowItems[j], details);
        if (rowItems[j].message) {
          img.attr('title', rowItems[j].message);
        }
      }
    }
    categoriesTable.append(row);
  }
}

function selectTutorialItems(items, expected) {
  var itemsByExpected = _.groupBy(items, function(x) { return x.expected; });
  var tutorial_items = [];
  for (var i = 0; i < expected.length; i++) {
    var items = itemsByExpected[expected[i]];
    if (items && items.length) {
      var si = Math.floor(Math.random()*items.length);
      tutorial_items.push(items[si]);
      items.splice(si, 1);
    }
  }
  return tutorial_items;
}


// Create UI for selecting category
// Returns array of ui elements corresponding to each category
function createSelectCategoryUI(parentElem, categories, keyCodeMap, disabled) {
  // Create UI for one item
  var uielems = [];
  for (var i = 0; i < categories.length; i++) {
    var category = categories[i];
    var inputElem = $("<input/>").attr("type", "radio").attr("name", category.label)
              .attr("value", category.value).prop("hidden", true).prop("disabled", disabled);
    var shortcut = (category.shortcut != null)? ' (' + category.shortcut + ')': '';
    var labelElem = $("<label></label>").attr("class", "btn btn-default")
          .text(category.label + shortcut) 
          .attr("title", category.description).append(inputElem);
    var categoryButton = $("<div></div>").attr("class", "radio-inline").append(labelElem).prop("disabled", disabled);
    if (category.addRow) {
      parentElem.append("<br/>");
    }
    parentElem.append(categoryButton);
    uielems.push({ input: inputElem, button: categoryButton });
    if (keyCodeMap != null && category.shortcut != null) {
      keyCodeMap[category.shortcut.toUpperCase().charCodeAt(0)] = i;
    }
  }
  return uielems;
}

function createItemImage(item, parent) {
  var imgElem = $('<img>').attr('src', item.image1)
    .data('alt-src', item.image2)
    .addClass('withHover')
    .addClass('bordered')
    .appendTo(parent);
  imgElem.hover(swap_image_source, swap_image_source);
  return imgElem;  
}

//
// Create symmetry task
// param categories {label: string, value: string, description: string, shortcut: string}
// param input {id: string, image1: string, image2: string, expected: string}
// param submitHookup {setupSubmit: function(), setOutput: function(output)}
// param useKeyshortcuts {boolean}
function SymmetryTask(categories, input, submitHookup, useKeyshortcuts) {
  this.input = input;
  // Annotations (parallel to input)
  this.annotations = [];
  for (var i = 0; i < this.input.length; i++) {
    this.annotations.push(null);
  }
  // Hookup for submit button
  this.submitHookup = submitHookup;
  // Whether to use keyboard shortcuts
  this.useKeyshortcuts = useKeyshortcuts;
  this.keyCodeMap = this.useKeyshortcuts? {} : null;

  this.categories = categories;
  this.valid_values = categories.map(function(x) { return x.value; });

  // Some variables to track state of the HIT.
  this.idx = 0;
  this.enabled = false;

  // UI elements
  this.numEntriesElem = $('#numEntries');
  this.imageContainer = $('#image-container');
  this.mainCategorySelectElem = $('#choices');
  this.counterTopElem = $('#counter .counter-top');
  this.counterBottomElem = $('#counter .counter-bottom');
  this.prevButton = $('#prev-btn');
  this.nextButton = $('#next-btn');
  this.submitButton = $('#submit-btn');
  // Track UI elements (corresponding to the categories)
  this.selectCategoryUIElems = [];
}

SymmetryTask.prototype.init = function(isPreview) {
  this.numEntriesElem.text(this.input.length);
  this.selectCategoryUIElems = createSelectCategoryUI(this.mainCategorySelectElem, this.categories, this.keyCodeMap, true);

  // Enable the UI if the HIT is not in preview mode.
  if (!isPreview) {
    this.enableHit();
  }

  // Set up the annotations.
  this.update();
};

SymmetryTask.prototype.showAnnotation = function(value) {
  this.mainCategorySelectElem.find('label').removeClass('active');
  this.mainCategorySelectElem.find('label input').prop('checked', false);
  if (value != null) {
    var inputElem = this.mainCategorySelectElem.find('label input[value=' + value + ']');
    inputElem.prop('checked', true);
    inputElem.parent().addClass('active');
  }
};

SymmetryTask.prototype.update = function() {
  // Set up the image
  var idx = this.idx;

  // Set up the image
  this.imageContainer.empty();
  createItemImage(this.input[idx], this.imageContainer);
  this.showAnnotation(this.annotations[idx]);
  this.updateSelected();

  // Refresh the counter
  this.counterTopElem.text(idx + 1);
  this.counterBottomElem.text(this.input.length);

  // If the UI is enabled, enable or disable the buttons depending on
  // the index.
  if (this.enabled) {
    var prev_btn = this.prevButton;
    var next_btn = this.nextButton;
    prev_btn.prop('disabled', true);
    next_btn.prop('disabled', true);
    if (idx > 0) {
      prev_btn.prop('disabled', false);
    }
    if (this.allowNext()) {
      next_btn.prop('disabled', false);
    }
  }
};

SymmetryTask.prototype.allowNext = function() {
  return this.idx < this.input.length - 1;
}

// Save current annotation
SymmetryTask.prototype.saveCurrent = function() {
  this.annotations[this.idx] = this.getSelected();
};

SymmetryTask.prototype.getSelected = function() {
  return this.mainCategorySelectElem.find('input:radio:checked').val();
};

SymmetryTask.prototype.updateSelected = function(input) {
  this.saveCurrent();
};

// Update the index, and save the text in the text area.
SymmetryTask.prototype.setIdx = function(new_idx) {
  if (new_idx < 0 || new_idx >= this.input.length) return;
  this.saveCurrent();

  this.idx = new_idx;
  this.update();
};

SymmetryTask.prototype.checkAnnotations = function() {
  var valid_values = this.valid_values;
  var invalidIndex = _.findIndex(this.annotations, function(d) { return d == null || valid_values.indexOf(d) < 0; });
  if (invalidIndex >= 0) {
    alert('Not all images annotated. Correct entry #' + (invalidIndex+1) + ' before submitting.');
    return false;
  }
  return true;
};

SymmetryTask.prototype.getOutput = function() {
  var output = _.map(_.zip(this.input, this.annotations), function(x) {
    return {'id': x[0].id, 'annotation': x[1]};
  });
  return output;
};

SymmetryTask.prototype.setOutput = function(output) {
  if (output.length !== this.input.length) {
    throw "Output does not match input";
  }
  var annotations = [];
  for (var i = 0; i < output.length; i++) {
    if (this.input[i].id !== output[i].id) {
      throw "Output id does not match input id for element " + i; 
    }
    annotations[i] = output[i].annotation;
  }
  this.annotations = annotations;
}

SymmetryTask.prototype.hookupKeys = function() {
  if (this.useKeyshortcuts) {
    var that = this;
    window.addEventListener("keyup", function (event) {
      var tagName = (event.target || event.srcElement).tagName;
      var isInputText = (tagName === 'INPUT' || tagName === 'SELECT' || tagName === 'TEXTAREA');
      if (!isInputText && that.enabled) {
        if (event.keyCode === 37 /* left arrow */) {
          that.prevButton.click();
        } else if (event.keyCode == 39 /* right arrow */) {
          that.nextButton.click();
        }
        var index = that.keyCodeMap[event.keyCode];
        if (index != null) {
          that.selectCategoryUIElems[index].input.click();
          that.selectCategoryUIElems[index].button.click();
          return false;
        }
      }
    });
  }
};

SymmetryTask.prototype.next = function() {
  this.setIdx(this.idx + 1);
};

SymmetryTask.prototype.prev = function() {
  this.setIdx(this.idx - 1);
};

SymmetryTask.prototype.enableHit = function() {
  this.enabled = true;

  // Enable components
  var that = this;
  // this.nextButton.click(function() { that.setIdx(that.idx + 1); });
  // this.prevButton.click(function() { that.setIdx(that.idx - 1); });
  this.nextButton.click(function() { that.next(); });
  this.prevButton.click(function() { that.prev(); });
  this.submitButton.prop('disabled', false);
  for (var i = 0; i < this.selectCategoryUIElems.length; i++) {
    var inputElem = this.selectCategoryUIElems[i];
    inputElem.input.prop('disabled', false);
    inputElem.button.prop('disabled', false);
    inputElem.button.click(function(event) {
      that.mainCategorySelectElem.find('label').removeClass('active');
      that.mainCategorySelectElem.find('label input').prop('checked', false);
      $(this).find('label').addClass('active');
      $(this).find('label input').prop('checked', true);
      that.updateSelected();
    });
  }
  this.hookupKeys();

  // Set up submit handler.
  this.submitHookup.setupSubmit();
  this.submitButton.click(function() {
    that.saveCurrent();
    var annotationsOkay = that.checkAnnotations();
    if (!annotationsOkay) { return false; }
    var output = that.getOutput();
    that.submitHookup.setOutput(output);
  });
};

function SymmetryTaskTutorial(categories, input) {
  // Tutorial for Symmetry Task
  var that = this;
  SymmetryTask.call(this, categories, input, { setupSubmit: function() {}, setOutput: function(output) { that.setOutput(output); } }, false);
  this.numEntriesElem = $('#tutorial-numEntries');
  this.imageContainer = $('#tutorial-image-container');
  this.mainCategorySelectElem = $('#tutorial-choices');
  this.counterTopElem = $('#tutorial-counter .counter-top');
  this.counterBottomElem = $('#tutorial-counter .counter-bottom');
  this.prevButton = $('#tutorial-prev-btn');
  this.nextButton = $('#tutorial-next-btn');
  this.submitButton = $('#tutorial-submit-btn');
  this.alertElem = $('#tutorial-message');
  this.mainElem = $('#tutorial-main');
  this.tutorialModal = $('#tutorial-modal');
  this.okay = false;
}

SymmetryTaskTutorial.prototype = Object.create(SymmetryTask.prototype);
SymmetryTaskTutorial.prototype.constructor = SymmetryTaskTutorial;

SymmetryTaskTutorial.prototype.showMessage = function(message, style, div) {
  div = div || $('<div></div>');
  div.attr('class', 'alert ' + style);
  div.empty();
  if (Array.isArray(message)) {
    for (var i = 0; i < message.length; i++) {
      if (i > 0) {
        div.append('<br/>');
      }
      div.append(message[i]);
    }
  } else {
    div.text(message);
  }
  return div;
};

SymmetryTaskTutorial.prototype.updateSelected = function() {
  SymmetryTask.prototype.updateSelected.call(this);
  var result = this.checkExpected(this.input[this.idx], this.annotations[this.idx]);
  //console.log('check', this.input[this.idx], this.annotations[this.idx], result);
  if (result.error) {
    this.showMessage(result.error, 'alert-danger', this.alertElem);
    this.okay = false;
  } else {
    this.showMessage(result.message, 'alert-success', this.alertElem);
    this.okay = true;
  }
  this.nextButton.prop('disabled', !this.okay);
};


SymmetryTaskTutorial.prototype.allowNext = function() {
  return this.okay;
};

SymmetryTaskTutorial.prototype.setIdx = function(new_idx) {
  if (new_idx >= this.input.length) {
    this.mainElem.text("Thank you for doing this tutorial.  Please press 'Start' to continue with the main task");
    this.submitButton.show();
    this.tutorialModal.find('.close').show();
  } else {
    SymmetryTask.prototype.setIdx.call(this, new_idx);
  }
};

SymmetryTaskTutorial.prototype.checkExpected = function(item, value) {
  var category = this.categories.find(function(x) { return x.value === item.expected; }) || {};
  if (value == null) {
    return { error: ["Move you cursor in and out of the image box to see the two images.", "Please indicate if you can detect any differences between the images"] };
  } else if (value === item.expected) {
    var message = 'Yes, that is correct! ' + (item.message || category.tutorial_message || '');
    return { message: message };
  } else {
    var error = item.error || "Hmm, that doesn't seem quite right";
    return { error: error };
  }
};


