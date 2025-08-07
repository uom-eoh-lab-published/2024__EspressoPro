#!/bin/bash
# -*- coding: utf-8 -*-
"""
Test runner script for EspressoPro
"""

echo "🧪 EspressoPro Test Suite"
echo "=========================="

# Change to the project directory
cd "$(dirname "$0")"

echo "📂 Current directory: $(pwd)"

# Run the main test runner
echo ""
echo "🔧 Running core functionality tests..."
python tests/run_tests.py

# Run individual tests
echo ""
echo "🔧 Running simple path detection test..."
python tests/test_simple.py

echo ""
echo "🔧 Running automatic paths test (if pytest available)..."
if command -v pytest &> /dev/null; then
    echo "✅ pytest found, running pytest..."
    pytest tests/ -v --tb=short
else
    echo "⚠️  pytest not found, running manual test..."
    python tests/test_automatic_paths.py
fi

echo ""
echo "🎉 Test suite completed!"
