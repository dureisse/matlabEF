#!/bin/csh

foreach i (`ls test*.m`)
  echo $i
  octave-cli $i
end
