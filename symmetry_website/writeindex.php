<?php
header("Cache-Control: no-store, no-cache, must-revalidate, max-age=0");
header("Cache-Control: post-check=0, pre-check=0", false);
header("Pragma: no-cache");


$doc = new DOMDocument('1.0');
$doc->formatOutput = true;

$root = $doc->createElement("html");
$root = $doc->appendChild($root);

$head = $doc->createElement("head");
$head = $root->appendChild($head);

$meta = $doc->createElement("meta");
$meta = $head->appendChild($meta);
$meta->setAttribute("http-equiv", "Pragma");
$meta->setAttribute("content", "no-cache");

$meta = $doc->createElement("meta");
$meta = $head->appendChild($meta);
$meta->setAttribute("http-equiv", "Cache-Control");
$meta->setAttribute("content", "no-cache");

$meta = $doc->createElement("meta");
$meta = $head->appendChild($meta);
$meta->setAttribute("http-equiv", "Expires");
$meta->setAttribute("content", "0");

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
$script->setAttribute("src", "indexcallbacks.js");
$script = $head->appendChild($script);

$body = $doc->createElement("body");
$body = $doc->appendChild($body);

$root = "/orions4-zfs/projects/mengyuan/julia_code/symmetry/mitsuba_render/results/";
$files = scandir($root);


foreach ($files as $file) {
    if (strpos($file, ".txt")) {
        $synset = explode(".", $file);
        $synset = $synset[0];
        $lines = count(file($root.$file));
        print_r($lines);
        $header = $doc->createElement("h1", $synset);
        $header = $body->appendChild($header);

        $table = $doc->createElement("table");
        $table->setAttribute("id", "T".$synset);
        $table = $body->appendChild($table);

        $numCell = ceil($lines / 45);
        for ($i = 1; $i <= ceil($numCell/10); $i++){
            $row = $doc->createElement("tr");
            $row = $table->appendChild($row);
            for ($j = 1; $j <= 10; $j++) {
                $index = ($i - 1) * 10 + $j;
                if ($index > $numCell) continue;
                $column = $doc->createElement("td");
                $column = $row->appendChild($column);
                $start = ($index - 1) * 45;
                $ending = $index * 45;
                if ($ending > $lines) $ending = $lines;
                $button = $doc->createElement("button", sprintf("%d -- %d", $start, $ending));
                $button->setAttribute("type", "button");
                $button->setAttribute("style", "width:250px; font-size:24px; background:'';");
                $button->setAttribute("data-synset", $synset);
                $button->setAttribute("data-start", $start);
                $button->setAttribute("data-end", $ending);
                $button = $column->appendChild($button);
            }
        }
    }
}


echo $doc->saveHTMLFile("index.html");

echo "done";

?>
