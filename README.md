# VOL<GO> Options Terminal

Static Bloomberg-style front end for the options pricing notebook at `C:\Users\c_sid\Downloads\Options_Pricing_Model (1).ipynb`.

## Included

- Terminal-inspired UI with command bar, pricing ticket, analytics ribbon, diagnostics, and activity feed
- Browser-side pricing engine modeled after the notebook's Black-Scholes, binomial tree, Monte Carlo, Greeks, and strategy logic
- Interactive payoff chart, synthetic volatility smile, and implied volatility heatmap with no build step

## Run

Open `C:\Users\c_sid\OneDrive\Documents\New project\index.html` in a browser.

## Current limitations

- This build is fully static because Node and Python are not available in the current workspace runtime.
- Live `yfinance` market data from the notebook is not connected yet, so the UI uses manual market inputs and synthetic vol surfaces.
- The Heston view is implemented as a fast smile-aware overlay rather than full browser-side calibration.
