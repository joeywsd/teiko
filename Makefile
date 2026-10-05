PYTHON ?= python3

.PHONY: setup pipeline dashboard

setup:
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) pipeline.py

dashboard:
	$(PYTHON) -m streamlit run dashboard.py
