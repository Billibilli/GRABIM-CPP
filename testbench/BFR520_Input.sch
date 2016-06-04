<Qucs Schematic 0.0.19>
<Properties>
  <View=-212,-70,1043,888,0.909091,109,251>
  <Grid=10,10,1>
  <DataSet=BFR520_Input.dat>
  <DataDisplay=BFR520_Input.dpl>
  <OpenDisplay=0>
  <Script=BFR520_Input.m>
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
  <Eqn Eqn3 1 550 50 -31 19 0 0 "S11dB=dB(S[1,1])" 1 "S21dB=dB(S[2,1])" 1 "S22dB=dB(S[2,2])" 1 "yes" 0>
  <GND * 1 180 410 0 0 0 0>
  <GND * 1 320 480 0 0 0 0>
  <.SP SP1 1 200 50 0 82 0 0 "lin" 1 "0.1 GHz" 1 "3 GHz" 1 "100" 1 "no" 0 "1" 0 "2" 0 "no" 0 "no" 0>
  <GND * 1 10 480 0 0 0 0>
  <R R1 1 320 430 15 -26 0 1 "50 Ohm" 1 "26.85" 0 "0.0" 0 "0.0" 0 "26.85" 0 "US" 0>
  <SPfile X1 1 180 380 -26 -67 0 0 "/home/andresmmera/git/GRABIM-CPP/testbench/Spar_BFR540_commonemitter8V30mA.s2p" 1 "rectangular" 0 "linear" 0 "open" 0 "2" 0>
  <Pac P1 1 10 430 18 -26 0 1 "1" 1 "75 Ohm" 1 "-5 dBm" 0 "1 GHz" 0 "26.85" 0>
</Components>
<Wires>
  <210 380 320 380 "" 0 0 0 "">
  <320 380 320 400 "" 0 0 0 "">
  <320 460 320 480 "" 0 0 0 "">
  <10 380 10 400 "" 0 0 0 "">
  <10 460 10 480 "" 0 0 0 "">
  <10 380 150 380 "" 0 0 0 "">
</Wires>
<Diagrams>
  <Rect 514 527 418 284 3 #c0c0c0 1 00 1 0 5e+08 3e+09 0 -20 5 30 1 -1 0.5 1 315 0 225 "" "" "">
	<"S21dB" #0000ff 0 3 0 0 0>
	<"S11dB" #ff0000 0 3 0 0 0>
	<"S22dB" #ff00ff 0 3 0 0 0>
  </Rect>
  <Smith 230 800 200 200 3 #c0c0c0 1 00 1 0 1 1 1 0 4 1 1 0 1 1 315 0 225 "" "" "">
	<"S[1,1]" #0000ff 0 3 0 0 0>
  </Smith>
</Diagrams>
<Paintings>
</Paintings>
