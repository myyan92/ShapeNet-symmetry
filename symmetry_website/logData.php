<?php

$synset = $_REQUEST["category"];
$startIdx = $_REQUEST["start"];
$endIdx = $_REQUEST["end"];
$annotation = $_REQUEST["annotation"];

$filepattern="results/%s_s%d_e%d.txt";
$filename = sprintf($filepattern, $synset, $startIdx, $endIdx);

$file = fopen($filename, 'w');
$linepattern = "%d\n";
foreach ($annotation as $value){
    $line = sprintf($linepattern, $value);
    fwrite($file, $line);
}
fclose($file);

// begin writing new webpage 
$doc = new DOMDocument('1.0');
$doc->formatOutput = true;

$root = $doc->createElement("html");
$root = $doc->appendChild($root);

$head = $doc->createElement("head");
$head = $root->appendChild($head);

$title = $doc->createElement("title");
$title = $head->appendChild($title);

$text = $doc->createTextNode("Thank you");
$text = $title->appendChild($text);

$body = $doc->createElement("body");
$body = $doc->appendChild($body);

$header = $doc->createElement("h1", "Thank you!");
$header = $body->appendChild($header);
$header->setAttribute("align", "center");

$content = $doc->createElement("p", "Your annotations are recorded");
$content->setAttribute("align", "center");
$content = $body->appendChild($content);
$log = $doc->createElement("p", $filename);
$log->setAttribute("align", "center");
$log = $body->appendChild($log);

$filename = "thankyou.html";
if (file_exists($filename))
	unlink($filename);
$doc->saveHTMLFile($filename);

echo "thankyou.html";

?>

