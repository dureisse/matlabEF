#!/bin/csh

foreach i (`ls test*.dgibi`)
  echo $i
  castem21 $i
end
foreach i (`ls test*.m`)
  echo $i
  octave-cli $i
end
