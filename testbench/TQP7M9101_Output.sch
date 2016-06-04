<Qucs Schematic 0.0.19>
<Properties>
  <View=0,-280,982,638,1,0,0>
  <Grid=10,10,1>
  <DataSet=TQP7M9101_Output.dat>
  <DataDisplay=TQP7M9101_Output.dpl>
  <OpenDisplay=0>
  <Script=TQP7M9101_Output.m>
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
  <GND * 1 190 160 0 0 0 0>
  <GND * 1 330 230 0 0 0 0>
  <.SP SP1 1 210 -200 0 82 0 0 "lin" 1 "0.1 GHz" 1 "3 GHz" 1 "100" 1 "no" 0 "1" 0 "2" 0 "no" 0 "no" 0>
  <GND * 1 20 230 0 0 0 0>
  <SPfile X1 1 190 130 -26 -67 0 0 "/home/andresmmera/git/GRABIM-CPP/testbench/TQP7M9101.s2p" 1 "rectangular" 0 "linear" 0 "open" 0 "2" 0>
  <Pac P1 1 330 180 18 -26 0 1 "1" 1 "50 Ohm" 1 "-5 dBm" 0 "1 GHz" 0 "26.85" 0>
  <R R1 1 20 180 15 -26 0 1 "50 Ohm" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <Eqn Eqn1 1 560 -200 -31 19 0 0 "S11dB=dB(S[1,1])" 1 "yes" 0>
</Components>
<Wires>
  <330 210 330 230 "" 0 0 0 "">
  <20 210 20 230 "" 0 0 0 "">
  <330 130 330 150 "" 0 0 0 "">
  <220 130 330 130 "" 0 0 0 "">
  <20 130 20 150 "" 0 0 0 "">
  <20 130 160 130 "" 0 0 0 "">
</Wires>
<Diagrams>
  <Rect 524 277 418 284 3 #c0c0c0 1 00 1 0 5e+08 3e+09 0 -20 5 30 1 -1 0.5 1 315 0 225 "" "" "">
	<"S21dB" #0000ff 0 3 0 0 0>
	<"S11dB" #ff0000 0 3 0 0 0>
	<"S22dB" #ff00ff 0 3 0 0 0>
  </Rect>
  <Smith 240 550 200 200 3 #c0c0c0 1 00 1 0 1 1 1 0 4 1 1 0 1 1 315 0 225 "" "" "">
	<"S[1,1]" #0000ff 0 3 0 0 0>
  </Smith>
</Diagrams>
<Paintings>
</Paintings>
