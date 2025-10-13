#!/usr/bin/env bash
echo "🧹 Cleaning Nextflow workspace..."
nextflow clean -f 
rm -rf output/work .nextflow .nextflow.log*
echo "✅ All cleaned!"