#!/usr/bin/env perl

# Date: 2013-05-02
# Author: Jeff Cullis

use strict;
use warnings;
use LWP::Simple;
use Text::Iconv;
use Spreadsheet::XLSX;
use Spreadsheet::ParseExcel; # for .xls
use Getopt::Long;
use utf8;

binmode(STDOUT, ":utf8"); # For e.g. degree symbols in input to render.

sub setupOptions {
    my  $url = 'http://biodiversity.agr.gc.ca/svn/sequencing/454_sample_summary.xls';
    my $options = {};
    GetOptions($options,
        'url|u=s',
        'excel_filename|x=s',
        'tab_filename|t=s',
        );
    $options->{url} ||= $url;
    if (!defined $options->{excel_filename} and $options->{url} =~ /\/([^\/]+)$/) { 
        $options->{excel_filename} = $1; 
    }
    if (!defined $options->{tab_filename} and $options->{excel_filename} =~ /^(.*)\.xlsx?$/) {
        $options->{tab_filename} = $1 . ".tab";
    }
    return $options;
}

sub downloadFile {
    my $url = shift;
    my $xlsx_filename = shift;

    getstore($url, $xlsx_filename);
}

sub parseXLSX {
    my $xlsx_filename = shift;
    my $tab_filename = shift;
        
    my $converter = Text::Iconv ->new('utf-8', 'windows-1251');
    my $excel = Spreadsheet::XLSX->new($xlsx_filename, $converter);
    open(FTAB, '>', $tab_filename) or die "Error: could not open file ${tab_filename} for writing.\n";
    binmode(FTAB, ":utf8");
    #foreach my $sheet (@{$excel->{Worksheet}}) {
    my $sheet = ${$excel->{Worksheet}}[0]; {
        $sheet->{MaxRow} ||= $sheet->{MinRow};
        for my $row ($sheet->{MinRow} .. $sheet->{MaxRow}) {
            $sheet->{MaxCol} ||= $sheet->{MinRow};
            for my $col ($sheet->{MinCol} .. $sheet->{MaxCol}) {
                my $cell = $sheet->{Cells}[$row][$col];
                my $val = ($cell ? $cell->{Val} : "");
                print FTAB $val . "\t";
            }
            print FTAB "\n";
        }
    }
    close FTAB;
}

sub getCellValue
{
    my $sheet = shift;
    my $row = shift;
    my $col = shift;
    my $cell;
    my $cv;
    $cell = $sheet->{Cells}[$row][$col];
    $cv = ( $cell ? $cell->{Val} : "" );
    return $cv;
}

sub rowIsEmpty
{
    my $sheet = shift;
    my $row = shift;
    my $col = $sheet->{MaxCol};
    my $cv;
    my $is_empty = 1;
    while ($col >= 0) {
        $cv = getCellValue($sheet, $row, $col);
        $col--;
        $is_empty = 0 if $cv =~ /\S/;
    }
    return $is_empty;
}

# Find last non-whitespace-only row, reporting excess whitespace.
sub getLastRow
{
    my $sheet = shift;
    my $row = $sheet->{MaxRow};
    while( rowIsEmpty($sheet,$row) ) { $row-- };
    return $row;
}

sub parseXLSWorksheet
{
    my $FOUT = shift;    
    my $worksheet = shift;

    my ( $row_min, $row_max ) = $worksheet->row_range();
    my ( $col_min, $col_max ) = $worksheet->col_range();
    $row_max = getLastRow($worksheet);

    for my $row ( $row_min .. $row_max ) {
           for my $col ( $col_min .. $col_max ) {
               my $cell = $worksheet->get_cell( $row, $col );
               my $val = ($cell ? $cell->value() : "");
            print $FOUT $val . "\t";
           }
           print $FOUT "\n";
    }
}    

sub parseFirstXLSWorksheet
{
    my $FOUT = shift;
    my $workbook = shift;    
    my @sheets = $workbook->worksheets();
    my $worksheet = $sheets[0];
    parseXLSWorksheet($FOUT, $worksheet);
}

sub parseAllXLSWorksheets
{
    my $FOUT = shift;
    my $workbook = shift;        
    for my $worksheet ( $workbook->worksheets() ) {
        parseXLSWorksheet($FOUT, $worksheet);
    }
}

sub parseXLS {
    my $excel_filename = shift;
    my $tab_filename = shift;
    
    my $parser   = Spreadsheet::ParseExcel->new();
    my $workbook = $parser->Parse($excel_filename);
    open(my $FTAB, '>', $tab_filename) or die "Error: could not open file ${tab_filename} for writing.\n";
    binmode($FTAB, ":utf8");
    # parseAllXLSWorksheets($FTAB, $workbook); # Note: outputs all sheets into one file - may not be desirable.
    parseFirstXLSWorksheet($FTAB, $workbook);
    close FTAB;
}

sub parseExcel {
    my $excel_filename = shift;
    my $tab_filename = shift;
    if ($excel_filename =~ /\.xls$/) {
        parseXLS ($excel_filename, $tab_filename);
    } elsif ($excel_filename =~ /\.xlsx$/) {
        parseXLSX ($excel_filename, $tab_filename);
    } else {
        die "Error: Excel input file ${excel_filename} does not have .xls or .xlsx extension\n";
    }
}

my $opts = setupOptions;
downloadFile($opts->{url}, $opts->{excel_filename});
parseExcel ($opts->{excel_filename}, $opts->{tab_filename});






