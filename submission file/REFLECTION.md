# Reflection

## why the app works
This app is based on the Michaelis-Menten equation, v = Vmax·[S] / (Km + [S]). Vmax is the maximum velocity an enzyme can achieve when all of its active sites are binded with substrate, which means it is saturated. While Km is the substrate concentration needed for an enzyme to reach half of its maximum reaction velocity(Vmax), it shows us how an enzyme interacts with its substrate
An inhibitor can change both Vmax and Km values by different ways of binding.
A competitive inhibitor binds the active site, so it competes with the substrate.When the substrate concentration is very high, the substrate wins over the inhibitor, thus the Vmax is unchanged, but more substrate is needed to reach the half Vmax, so there is an increase in Km value.
A non-competitive inhibitor binds to an allosteric site on an enzyme, which can slow down the enzyme's activity and cause reduction of the final product, regardless how much substrate there is.Thus Vmax decrease while Km unchanged.
A mixed inhibitor does both.
The app works because it measures Vmax and Km for each condition and compares them with the control using these rules

## Why the app works: The methodology
As the app fits the curve directly with R's `nls()`, instead of using a lineweaver-Burk plot.Taking reciprocal values (1/v and 1/[S]), which may ampilfy small errors at low substrate concentration into large errors and distort the result. Fitting the original curve directly can avoid this problem and gives a reliable standard error for each value.
i think the app works because i have tested it with data that i already knew the answer.For the competitive test, the app found Vmax unchanged(100.17 V.S. 100.8) and Km increased(2.008 V.S. 5.08).For the second test, the app correctly labelled non-competitive and mixed inhibitor. Also, when i removed a value(using the test_missing.csv), the app will appear a messang that warnd me, and fitted the remaining 6 points instead of crashing

## Where it could fail
The app assumes the first condition in the file is the control, which means a differently ordered dataset may produce an invalid comparison.
Also, this app cannot use to measure enzyme that deviate from standard Michaelis–Menten kinetics, such as cooperative enzymes or those exhibiting substrate inhibition.
Besides, i haven't tested the app with uncompetitive inhibition or cases where the nls() fit fails to converge.

## Working with the AI agent
The biggest lesson i have learnt is about the prompts used to communicate with the AI should not be ambigous.For example, i asked for rounding 'to 2 significant figures' and the AI agent applied this to everything, giving "1.0e2 ± 0.96", which i. don't want.I realised the problem is my instruction was a bit ambigous, thus, i wrote the instruction more precisely and give examples in the second attempt,it finally worked.
Also, something surprising me, the example i gave said the value is 100.20, but the real estimate was 100.17. The AI agent followed the rule instead of my wrong example, thus, this app was right and i was wrong. This shows me that checking has to be go both ways, I cannot just trust the AI, meanwhile, i cannot assume my own expectations are all correct. The only reliable test is to check the real numbers

## Conclusion
I will trust this app to give me a quick and first look at inhibition data. However, i should not treat this as a final answer for the real lab results without any further checking, such as a proper statistical comparison.