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
    parentElem.append(categoryButton);
    uielems.push({ input: inputElem, button: categoryButton });
    if (keyCodeMap != null && category.shortcut != null) {
      keyCodeMap[category.shortcut.toUpperCase().charCodeAt(0)] = i;
    }
  }
  return uielems;
}

//
// Create symmetry task
// param categories {label: string, value: string, description: string, shortcut: string}
// param input {id: string, image1: string, image2: string, expected: string}
// param submitHookup {setupSubmit: function(), setOutput: function(output)}
// param useKeyshortcuts {boolean}
//
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
  this.counterTopElem = $('.counter-top');
  this.counterBottomElem = $('.counter-bottom');
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
    this.enable_hit();
  }

  this.preloadImages();

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

SymmetryTask.prototype.renderItem = function(item, annotation) {
  // Set up the image
  this.imageContainer.empty();
  var imgElem = $('<img>').attr('src', item.image1)
    .data('alt-src', item.image2)
    .addClass('withHover')
    .addClass('bordered')
    .appendTo(this.imageContainer);
  imgElem.hover(swap_image_source, swap_image_source);
  this.showAnnotation(annotation);
};

SymmetryTask.prototype.update = function() {
  // Set up the image
  var idx = this.idx;
  this.renderItem(this.input[idx], this.annotations[idx])

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
    if (idx < this.input.length - 1) {
      next_btn.prop('disabled', false);
    }
  }
};

SymmetryTask.prototype.preloadImages = function() {
  preloadImages(this.input.map(function(x) { return x.image1; }));
  preloadImages(this.input.map(function(x) { return x.image2; }));
}

// Save current annotation
SymmetryTask.prototype.save_current = function() {
  this.annotations[this.idx] = this.mainCategorySelectElem.find('input:radio:checked').val();
};

// Update the index, and save the text in the text area.
SymmetryTask.prototype.set_idx = function(new_idx) {
  if (new_idx < 0 || new_idx >= this.input.length) return;
  this.save_current();

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

SymmetryTask.prototype.hookupKeys = function() {
  if (this.useKeyshortcuts) {
    var that = this;
    window.addEventListener("keyup", function (event) {
      var tagName = (event.target || event.srcElement).tagName;
      var isInputText = (tagName === 'INPUT' || tagName === 'SELECT' || tagName === 'TEXTAREA');
      if (!isInputText) {
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

SymmetryTask.prototype.enable_hit = function() {
  this.enabled = true;

  // Enable components
  var that = this;
  this.nextButton.click(function() { that.set_idx(that.idx + 1); });
  this.prevButton.click(function() { that.set_idx(that.idx - 1); });
  this.submitButton.prop('disabled', false);
  for (var i = 0; i < this.selectCategoryUIElems.length; i++) {
    var inputElem = this.selectCategoryUIElems[i];
    inputElem.input.prop('disabled', false);
    inputElem.button.prop('disabled', false);
    inputElem.button.click(function(event) {
      that.mainCategorySelectElem.find('label').removeClass('active');
      $(this).find('label').addClass('active');
    });
  }
  this.hookupKeys();

  // Set up submit handler.
  this.submitHookup.setupSubmit();
  this.submitButton.click(function() {
    that.save_current();
    var annotationsOkay = that.checkAnnotations();
    if (!annotationsOkay) { return false; }
    var output = that.getOutput();
    that.submitHookup.setOutput(output);
  });
}
