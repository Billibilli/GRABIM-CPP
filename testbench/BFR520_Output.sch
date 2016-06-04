<Qucs Schematic 0.0.19>
<Properties>
  <View=0,-267,1062,625,1,0,0>
  <Grid=10,10,1>
  <DataSet=BFR520_Output.dat>
  <DataDisplay=BFR520_Output.dpl>
  <OpenDisplay=0>
  <Script=BFR520_Output.m>
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
  <GND * 1 270 190 0 0 0 0>
  <GND * 1 410 260 0 0 0 0>
  <.SP SP1 1 290 -170 0 82 0 0 "lin" 1 "0.1 GHz" 1 "3 GHz" 1 "100" 1 "no" 0 "1" 0 "2" 0 "no" 0 "no" 0>
  <GND * 1 100 260 0 0 0 0>
  <SPfile X1 1 270 160 -26 -67 0 0 "/home/andresmmera/git/GRABIM-CPP/testbench/Spar_BFR540_commonemitter8V30mA.s2p" 1 "rectangular" 0 "linear" 0 "open" 0 "2" 0>
  <Pac P1 1 410 210 18 -26 0 1 "1" 1 "50 Ohm" 1 "-5 dBm" 0 "1 GHz" 0 "26.85" 0>
  <R R1 1 100 210 15 -26 0 1 "75 Ohm" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <Eqn Eqn1 1 640 -170 -31 19 0 0 "S11dB=dB(S[1,1])" 1 "yes" 0>
</Components>
<Wires>
  <300 160 410 160 "" 0 0 0 "">
  <410 160 410 180 "" 0 0 0 "">
  <410 240 410 260 "" 0 0 0 "">
  <100 240 100 260 "" 0 0 0 "">
  <100 160 240 160 "" 0 0 0 "">
  <100 160 100 180 "" 0 0 0 "">
</Wires>
<Diagrams>
  <Rect 604 307 418 284 3 #c0c0c0 1 00 1 0 5e+08 3e+09 0 -20 5 30 1 -1 0.5 1 315 0 225 "" "" "">
	<"S21dB" #0000ff 0 3 0 0 0>
	<"S11dB" #ff0000 0 3 0 0 0>
	<"S22dB" #ff00ff 0 3 0 0 0>
  </Rect>
  <Smith 320 580 200 200 3 #c0c0c0 1 00 1 0 1 1 1 0 4 1 1 0 1 1 315 0 225 "" "" "">
	<"S[1,1]" #0000ff 0 3 0 0 0>
  </Smith>
</Diagrams>
<Paintings>
</Paintings>
