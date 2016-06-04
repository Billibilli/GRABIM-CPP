<Qucs Schematic 0.0.19>
<Properties>
  <View=0,-197,962,800,1,0,0>
  <Grid=10,10,1>
  <DataSet=BFG310_Output.dat>
  <DataDisplay=BFG310_Output.dpl>
  <OpenDisplay=0>
  <Script=BFG310_Output.m>
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
  <Eqn Eqn1 1 580 -180 -31 19 0 0 "S11dB=dB(S[1,1])" 1 "S21dB=dB(S[2,1])" 1 "S22dB=dB(S[2,2])" 1 "yes" 0>
  <GND * 1 210 180 0 0 0 0>
  <GND * 1 350 250 0 0 0 0>
  <.SP SP1 1 230 -180 0 82 0 0 "lin" 1 "0.1 GHz" 1 "3 GHz" 1 "100" 1 "no" 0 "1" 0 "2" 0 "no" 0 "no" 0>
  <GND * 1 40 250 0 0 0 0>
  <SPfile X1 1 210 150 -26 -67 0 0 "/home/andresmmera/git/GRABIM-CPP/testbench/Spar_BFG310WXR_3V5mA.s2p" 1 "rectangular" 0 "linear" 0 "open" 0 "2" 0>
  <R R1 1 40 200 15 -26 0 1 "50 Ohm" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <Pac P1 1 350 200 18 -26 0 1 "1" 1 "50 Ohm" 1 "-5 dBm" 0 "1 GHz" 0 "26.85" 0>
</Components>
<Wires>
  <240 150 350 150 "" 0 0 0 "">
  <350 150 350 170 "" 0 0 0 "">
  <350 230 350 250 "" 0 0 0 "">
  <40 150 40 170 "" 0 0 0 "">
  <40 230 40 250 "" 0 0 0 "">
  <40 150 180 150 "" 0 0 0 "">
</Wires>
<Diagrams>
  <Rect 544 297 418 284 3 #c0c0c0 1 00 1 0 5e+08 3e+09 0 -20 5 30 1 -1 0.5 1 315 0 225 "" "" "">
	<"S21dB" #0000ff 0 3 0 0 0>
	<"S11dB" #ff0000 0 3 0 0 0>
	<"S22dB" #ff00ff 0 3 0 0 0>
  </Rect>
  <Smith 260 570 200 200 3 #c0c0c0 1 00 1 0 1 1 1 0 4 1 1 0 1 1 315 0 225 "" "" "">
	<"S[1,1]" #0000ff 0 3 0 0 0>
  </Smith>
</Diagrams>
<Paintings>
</Paintings>
