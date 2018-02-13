<?php
header("Cache-Control: no-store, no-cache, must-revalidate, max-age=0");
header("Cache-Control: post-check=0, pre-check=0", false);
header("Pragma: no-cache");

$synset = $_REQUEST["category"];
$startIdx = $_REQUEST["start"];
$endIdx = $_REQUEST["end"];
$result_name = sprintf("results/%s_s%s_e%s.txt", $synset, $startIdx, $endIdx);
if (file_exists($result_name)) {
    echo "Done";
} else {
$modellist = file("/orions4-zfs/projects/mengyuan/julia_code/symmetry/mitsuba_render/results/".$synset.".txt", FILE_IGNORE_NEW_LINES | FILE_SKIP_EMPTY_LINES);
$imageroot = "media/".$synset."/";

$doc = new DOMDocument('1.0');
$doc->formatOutput = true;

$root = $doc->createElement("html");
$root = $doc->appendChild($root);

$head = $doc->createElement("head");
$head = $root->appendChild($head);

$title = $doc->createElement("title");
$title = $head->appendChild($title);

$text = $doc->createTextNode("symmetry detection verification");
$text = $title->appendChild($text);

$script = $doc->createElement("script");
$script->setAttribute("type", "text/javascript");
$script->setAttribute("src", "https://ajax.googleapis.com/ajax/libs/jquery/1.10.2/jquery.min.js");
$script = $head->appendChild($script);
$script = $doc->createElement("script");
$script->setAttribute("type", "text/javascript");
$script->setAttribute("src", "callbacks.js");
$script = $head->appendChild($script);

$body = $doc->createElement("body");
$body = $doc->appendChild($body);

$header = $doc->createElement("h1", "Symmetry results");
$header = $body->appendChild($header);
$header->setAttribute("align", "center");

$table = $doc->createElement("table");
$table = $body->appendChild($table);

for ($i = 1; $i <= 15; $i++){
  $row = $doc->createElement("tr");
  $row = $table->appendChild($row);
  for ($j = 1; $j <= 3; $j++) {
    $column = $doc->createElement("td");
    $column = $row->appendChild($column);
    $subTable = $doc->createElement("table");
    $subTable = $column->appendChild($subTable);
    $idx = ($i-1)*3+($j-1) + $startIdx;
    if ($idx >= $endIdx) continue;
    $imagenameBase = $imageroot.$modellist[$idx];
    $subRow = $doc->createElement("tr");
    $subRow = $subTable->appendChild($subRow);
    $subColumn = $doc->createElement("td");
    $subColumn = $subRow->appendChild($subColumn);
    $subdiv = $doc->createElement("div");   // need a div for positioning images
    $subdiv->setAttribute("style", "position:relative");
    $subdiv = $subColumn->appendChild($subdiv);
    $image = $doc->createElement("img");
    $image->setAttribute("class", "withHover");
    $image->setAttribute("src", $imagenameBase."_v0.0.png");
    $image->setAttribute("data-alt-src", $imagenameBase."_v0.1.png");
    $image = $subdiv->appendChild($image);
    $icon = $doc->createElement("img");
    $icon->setAttribute("src", "check.jpeg");
    $icon->setAttribute("style", "position:absolute; top:0px; left:0px;");
    $icon = $subdiv->appendChild($icon);
    $subColumn = $doc->createElement("td");
    $subColumn = $subRow->appendChild($subColumn);
    $subdiv = $doc->createElement("div");   // to make structure symmetric
    $subdiv->setAttribute("style", "position:relative");
    $subdiv = $subColumn->appendChild($subdiv);
    $image = $doc->createElement("img");
    $image->setAttribute("class", "withHover");
    $image->setAttribute("src", $imagenameBase."_v1.0.png");
    $image->setAttribute("data-alt-src", $imagenameBase."_v1.1.png");
    $image = $subdiv->appendChild($image);

  }
}

$button = $doc->createElement("button", "Submit");
$button->setAttribute("type", "button");
$button->setAttribute("style", "height:100px; width:300px;font-size:28px");
$button->setAttribute("data-synset", $synset);
$button->setAttribute("data-start", $startIdx);
$button->setAttribute("data-end", $endIdx);
$button = $body->appendChild($button);

$htmlname = sprintf("%s_s%s_e%s.html", $synset, $startIdx, $endIdx);
$doc->saveHTMLFile($htmlname);
echo $htmlname;
}
?>
