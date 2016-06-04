<Qucs Schematic 0.0.19>
<Properties>
  <View=0,3,975,630,1,0,0>
  <Grid=10,10,1>
  <DataSet=test2.dat>
  <DataDisplay=test2.dpl>
  <OpenDisplay=0>
  <Script=test2.m>
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
  <GND * 1 130 120 0 0 0 0>
  <L L1 1 280 90 -26 10 0 0 "2.15863e-08" 1 "" 0>
  <GND * 1 340 210 0 0 0 0>
  <L L2 1 340 180 -26 2 0 1 "1.95563e-08" 1 "" 0>
  <C C1 1 420 90 -26 17 0 0 "6.27691e-12" 1 "" 0 "neutral" 0>
  <GND * 1 480 210 0 0 0 0>
  <C C2 1 480 180 -26 17 0 1 "3.37685e-13" 1 "" 0 "neutral" 0>
  <C C3 1 560 90 -26 17 0 0 "7.98717e-12" 1 "" 0 "neutral" 0>
  <GND * 1 620 210 0 0 0 0>
  <L L3 1 620 180 -26 2 0 1 "1.4318e-07" 1 "" 0>
  <Pac P1 1 670 180 18 -26 0 1 "1" 1 "50 Ohm" 1 "0 dBm" 0 "1 GHz" 0 "26.85" 0>
  <GND * 1 670 210 0 0 0 0>
  <.SP SP1 1 89 413 0 82 0 0 "lin" 1 "0.1 GHz" 1 "3 GHz" 1 "100" 1 "no" 0 "1" 0 "2" 0 "no" 0 "no" 0>
  <Eqn Eqn1 1 329 423 -31 19 0 0 "S11dB=dB(S[1,1])" 1 "yes" 0>
  <SPfile X1 1 130 90 -26 -67 1 2 "/home/andresmmera/git/GRABIM-CPP/testbench//BFG310-output.s1p" 1 "rectangular" 0 "linear" 0 "open" 0 "1" 0>
</Components>
<Wires>
  <160 90 250 90 "" 0 0 0 "">
  <310 90 340 90 "" 0 0 0 "">
  <340 90 340 150 "" 0 0 0 "">
  <340 90 390 90 "" 0 0 0 "">
  <450 90 480 90 "" 0 0 0 "">
  <480 90 480 150 "" 0 0 0 "">
  <480 90 530 90 "" 0 0 0 "">
  <590 90 620 90 "" 0 0 0 "">
  <620 90 620 150 "" 0 0 0 "">
  <670 90 670 150 "" 0 0 0 "">
  <620 90 670 90 "" 0 0 0 "">
</Wires>
<Diagrams>
  <Rect 498 579 437 296 3 #c0c0c0 1 00 1 0 0.2 1 1 -0.1 0.5 1.1 1 -0.1 0.5 1.1 315 0 225 "" "" "">
	<"S11dB" #0000ff 0 3 0 0 0>
  </Rect>
</Diagrams>
<Paintings>
</Paintings>
