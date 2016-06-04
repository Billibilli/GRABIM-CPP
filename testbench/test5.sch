<Qucs Schematic 0.0.19>
<Properties>
  <View=-44,-47,1113,800,1,0,47>
  <Grid=10,10,1>
  <DataSet=test5.dat>
  <DataDisplay=test5.dpl>
  <OpenDisplay=0>
  <Script=test5.m>
  <RunScript=0>
  <showFrame=0>
  <FrameText0=Title>
  <FrameText1=Drawn By:>
  <FrameText2=Date:>
  <FrameText3=Revision:>
</Properties>
<Symbol>
</Symbol>
<Components>
  <.SP SP1 1 259 353 0 82 0 0 "lin" 1 "0.1 GHz" 1 "3 GHz" 1 "100" 1 "no" 0 "1" 0 "2" 0 "no" 0 "no" 0>
  <Eqn Eqn1 1 499 363 -31 19 0 0 "S11dB=dB(S[1,1])" 1 "yes" 0>
  <GND * 1 90 110 0 0 0 0>
  <GND * 1 180 200 0 0 0 0>
  <C C1 1 180 170 -26 17 0 1 "1.1504e-12" 1 "" 0 "neutral" 0>
  <C C2 1 260 80 -26 17 0 0 "2.97701e-12" 1 "" 0 "neutral" 0>
  <GND * 1 320 200 0 0 0 0>
  <L L1 1 320 170 -26 2 0 1 "7.61154e-09" 1 "" 0>
  <Pac P1 1 370 170 18 -26 0 1 "1" 1 "50 Ohm" 1 "0 dBm" 0 "1 GHz" 0 "26.85" 0>
  <GND * 1 370 200 0 0 0 0>
  <SPfile X1 1 90 80 -26 -67 1 2 "/home/andresmmera/git/GRABIM-CPP/testbench/TQP7M9101-output.s1p" 1 "rectangular" 0 "linear" 0 "open" 0 "1" 0>
</Components>
<Wires>
  <120 80 180 80 "" 0 0 0 "">
  <180 80 180 140 "" 0 0 0 "">
  <180 80 230 80 "" 0 0 0 "">
  <290 80 320 80 "" 0 0 0 "">
  <320 80 320 140 "" 0 0 0 "">
  <370 80 370 140 "" 0 0 0 "">
  <320 80 370 80 "" 0 0 0 "">
</Wires>
<Diagrams>
  <Rect 668 519 437 296 3 #c0c0c0 1 00 1 0 0.2 1 1 -0.1 0.5 1.1 1 -0.1 0.5 1.1 315 0 225 "" "" "">
	<"S11dB" #0000ff 0 3 0 0 0>
  </Rect>
</Diagrams>
<Paintings>
</Paintings>
