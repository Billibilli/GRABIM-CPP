<Qucs Schematic 0.0.19>
<Properties>
  <View=-134,-60,1055,800,0.826446,0,37>
  <Grid=10,10,1>
  <DataSet=test4.dat>
  <DataDisplay=test4.dpl>
  <OpenDisplay=0>
  <Script=test4.m>
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
  <.SP SP1 1 169 363 0 82 0 0 "lin" 1 "0.1 GHz" 1 "3 GHz" 1 "100" 1 "no" 0 "1" 0 "2" 0 "no" 0 "no" 0>
  <Eqn Eqn1 1 409 373 -31 19 0 0 "S11dB=dB(S[1,1])" 1 "yes" 0>
  <GND * 1 460 100 0 0 0 0>
  <Pac P1 1 -80 160 18 -26 0 1 "1" 1 "50 Ohm" 1 "0 dBm" 0 "1 GHz" 0 "26.85" 0>
  <GND * 1 -80 190 0 0 0 0>
  <L L1 1 40 70 -26 10 0 0 "4.94551e-09" 1 "" 0>
  <C C1 1 160 70 -26 17 0 0 "2.16507e-11" 1 "" 0 "neutral" 0>
  <GND * 1 220 190 0 0 0 0>
  <C C2 1 220 160 -26 17 0 1 "1.41097e-12" 1 "" 0 "neutral" 0>
  <GND * 1 240 190 0 0 0 0>
  <L L2 1 240 160 -26 2 0 1 "3.81277e-09" 1 "" 0>
  <C C3 1 320 70 -26 17 0 0 "8.01808e-13" 1 "" 0 "neutral" 0>
  <GND * 1 380 190 0 0 0 0>
  <L L3 1 380 160 -26 2 0 1 "9.38882e-09" 1 "" 0>
  <GND * 1 400 190 0 0 0 0>
  <C C4 1 400 160 -26 17 0 1 "2.94952e-12" 1 "" 0 "neutral" 0>
  <SPfile X1 1 460 70 -26 -67 0 0 "/home/andresmmera/git/GRABIM-CPP/testbench/TQP7M9101-input.s1p" 1 "rectangular" 0 "linear" 0 "open" 0 "1" 0>
</Components>
<Wires>
  <-80 70 -80 130 "" 0 0 0 "">
  <-80 70 10 70 "" 0 0 0 "">
  <70 70 130 70 "" 0 0 0 "">
  <190 70 220 70 "" 0 0 0 "">
  <220 70 220 130 "" 0 0 0 "">
  <220 70 240 70 "" 0 0 0 "">
  <240 70 240 130 "" 0 0 0 "">
  <240 70 290 70 "" 0 0 0 "">
  <350 70 380 70 "" 0 0 0 "">
  <380 70 380 130 "" 0 0 0 "">
  <380 70 400 70 "" 0 0 0 "">
  <400 70 400 130 "" 0 0 0 "">
  <400 70 430 70 "" 0 0 0 "">
</Wires>
<Diagrams>
  <Rect 578 529 437 296 3 #c0c0c0 1 00 1 0 0.2 1 1 -0.1 0.5 1.1 1 -0.1 0.5 1.1 315 0 225 "" "" "">
	<"S11dB" #0000ff 0 3 0 0 0>
  </Rect>
</Diagrams>
<Paintings>
</Paintings>
