<Qucs Schematic 0.0.19>
<Properties>
  <View=0,-7,1101,620,1,0,60>
  <Grid=10,10,1>
  <DataSet=test3.dat>
  <DataDisplay=test3.dpl>
  <OpenDisplay=0>
  <Script=test3.m>
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
  <Pac P1 1 60 170 18 -26 0 1 "1" 1 "75 Ohm" 1 "0 dBm" 0 "1 GHz" 0 "26.85" 0>
  <GND * 1 60 200 0 0 0 0>
  <L L1 1 180 80 -26 10 0 0 "7.47241e-09" 1 "" 0>
  <C C1 1 300 80 -26 17 0 0 "2.29278e-12" 1 "" 0 "neutral" 0>
  <GND * 1 360 200 0 0 0 0>
  <L L2 1 360 170 -26 2 0 1 "5.43933e-09" 1 "" 0>
  <C C2 1 440 80 -26 17 0 0 "3.95189e-12" 1 "" 0 "neutral" 0>
  <GND * 1 500 200 0 0 0 0>
  <L L3 1 500 170 -26 2 0 1 "1.82325e-09" 1 "" 0>
  <GND * 1 520 200 0 0 0 0>
  <C C3 1 520 170 -26 17 0 1 "9.3985e-12" 1 "" 0 "neutral" 0>
  <GND * 1 570 110 0 0 0 0>
  <SPfile X1 1 570 80 -26 -67 0 0 "/home/andresmmera/git/GRABIM-CPP/testbench/BFR520-input.s1p" 1 "rectangular" 0 "linear" 0 "open" 0 "1" 0>
  <.SP SP1 1 79 373 0 82 0 0 "lin" 1 "0.1 GHz" 1 "3 GHz" 1 "100" 1 "no" 0 "1" 0 "2" 0 "no" 0 "no" 0>
  <Eqn Eqn1 1 319 383 -31 19 0 0 "S11dB=dB(S[1,1])" 1 "yes" 0>
</Components>
<Wires>
  <60 80 60 140 "" 0 0 0 "">
  <60 80 150 80 "" 0 0 0 "">
  <210 80 270 80 "" 0 0 0 "">
  <330 80 360 80 "" 0 0 0 "">
  <360 80 360 140 "" 0 0 0 "">
  <360 80 410 80 "" 0 0 0 "">
  <470 80 500 80 "" 0 0 0 "">
  <500 80 500 140 "" 0 0 0 "">
  <500 80 520 80 "" 0 0 0 "">
  <520 80 520 140 "" 0 0 0 "">
  <520 80 540 80 "" 0 0 0 "">
</Wires>
<Diagrams>
  <Rect 488 539 437 296 3 #c0c0c0 1 00 1 0 0.2 1 1 -0.1 0.5 1.1 1 -0.1 0.5 1.1 315 0 225 "" "" "">
	<"S11dB" #0000ff 0 3 0 0 0>
  </Rect>
</Diagrams>
<Paintings>
</Paintings>
