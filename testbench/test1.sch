<Qucs Schematic 0.0.19>
<Properties>
  <View=0,13,1122,800,0.909091,0,42>
  <Grid=10,10,1>
  <DataSet=test1.dat>
  <DataDisplay=test1.dpl>
  <OpenDisplay=0>
  <Script=test1.m>
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
  <.SP SP1 1 199 413 0 82 0 0 "lin" 1 "0.1 GHz" 1 "3 GHz" 1 "100" 1 "no" 0 "1" 0 "2" 0 "no" 0 "no" 0>
  <Eqn Eqn1 1 439 423 -31 19 0 0 "S11dB=dB(S[1,1])" 1 "yes" 0>
  <Pac P1 1 80 210 18 -26 0 1 "1" 1 "50 Ohm" 1 "0 dBm" 0 "1 GHz" 0 "26.85" 0>
  <GND * 1 80 240 0 0 0 0>
  <GND * 1 140 240 0 0 0 0>
  <L L1 1 140 210 -26 2 0 1 "2.42845e-08" 1 "" 0>
  <GND * 1 160 240 0 0 0 0>
  <C C1 1 160 210 -26 17 0 1 "3.72211e-12" 1 "" 0 "neutral" 0>
  <L L2 1 240 120 -26 10 0 0 "9.02665e-09" 1 "" 0>
  <GND * 1 300 240 0 0 0 0>
  <C C2 1 300 210 -26 17 0 1 "5.77399e-12" 1 "" 0 "neutral" 0>
  <GND * 1 320 240 0 0 0 0>
  <L L3 1 320 210 -26 2 0 1 "1.45329e-08" 1 "" 0>
  <L L4 1 400 120 -26 10 0 0 "1.29449e-08" 1 "" 0>
  <GND * 1 490 150 0 0 0 0>
  <SPfile X3 1 490 120 -26 -67 0 0 "/home/andresmmera/git/GRABIM-CPP/testbench/BFG310-input.s1p" 1 "rectangular" 0 "linear" 0 "open" 0 "1" 0>
</Components>
<Wires>
  <80 120 80 180 "" 0 0 0 "">
  <80 120 140 120 "" 0 0 0 "">
  <140 120 140 180 "" 0 0 0 "">
  <140 120 160 120 "" 0 0 0 "">
  <160 120 160 180 "" 0 0 0 "">
  <160 120 210 120 "" 0 0 0 "">
  <270 120 300 120 "" 0 0 0 "">
  <300 120 300 180 "" 0 0 0 "">
  <300 120 320 120 "" 0 0 0 "">
  <320 120 320 180 "" 0 0 0 "">
  <320 120 370 120 "" 0 0 0 "">
  <430 120 460 120 "" 0 0 0 "">
</Wires>
<Diagrams>
  <Rect 608 579 437 296 3 #c0c0c0 1 00 1 0 0.2 1 1 -0.1 0.5 1.1 1 -0.1 0.5 1.1 315 0 225 "" "" "">
	<"S11dB" #0000ff 0 3 0 0 0>
  </Rect>
</Diagrams>
<Paintings>
</Paintings>
