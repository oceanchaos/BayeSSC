<?php
$fp=fopen("bayes.count","a");
if(!$fp) echo("Counter error")
else fputs($fp, "\n".date("m/d/y H:i:s"));
fclose($fp);
?>

