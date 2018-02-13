<?php
header("Cache-Control: no-store, no-cache, must-revalidate, max-age=0");
header("Cache-Control: post-check=0, pre-check=0", false);
header("Pragma: no-cache");

$synset = $_REQUEST["category"];
$startIdx = $_REQUEST["start"];
$endIdx = $_REQUEST["end"];
$result_name = sprintf("results/%s_s%s_e%s.txt", $synset, $startIdx, $endIdx);
if (file_exists($result_name)) {
    echo "exist";
} else {
    echo "None";
}
?>
