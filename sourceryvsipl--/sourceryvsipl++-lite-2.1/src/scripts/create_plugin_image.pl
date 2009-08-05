#! /usr/bin/perl
#########################################################################
# create_plugin_image.pl -- create binary plugin image from objdump	#
# Jules Bergmann, CodeSourcery, Inc					#
# Feb 12, 2009								#
#########################################################################

use strict;

sub create_plugin_image {
   my ($inf, $outf) = @_;

   open(OUT, "> $outf") and binmode OUT || die "Can't write $outf: $!\n";

   # 1. Find addresses of functions

   open(IN, $inf) || die "Can't read $inf: $!\n";

   my %func;
   while (<IN>) {
      if (/([0-9a-f]+) <(.+)>:/) {
         # printf("$2 - $1\n");
	 $func{$2} = hex($1);
	 }
      }
   close(IN);

   # printf("kernel: %08x\n", $func{"kernel"});
   # printf("input : %08x\n", $func{"input"});
   # printf("output: %08x\n", $func{"output"});

   # 2. Write out header block.

   print OUT pack 'N', $func{"kernel"};
   print OUT pack 'N', $func{"input"};
   print OUT pack 'N', $func{"output"};
   for (my $i=0; $i<32-3; $i+=1) {
     print OUT pack 'N*', 0;
     }

   # 3. Translate opcodes/data

   open(IN, $inf) || die "Can't read $inf: $!\n";

   my $section;
   my $size = 0;
   my $last_addr = hex("5000") - 4;

   my $conv = 0;

   while (<IN>) {
      if (/Disassembly of section \.(.+):/) {
         $section = $1;

	 if ($section =~ /^text/ ||
	     $section =~ /^data/ ||
	     $section =~ /^rodata/ ||
	     $section =~ /^ctors/) {
     	    $conv = 1;
	    }
	 else {
     	   $conv = 0;
	   }

	 # print "section: $section $conv\n";

	 next;
	 }
      next if $conv == 0;

      if (/^\s*([0-9a-f]+):\s*([0-9a-f][0-9a-f]) ([0-9a-f][0-9a-f]) ([0-9a-f][0-9a-f]) ([0-9a-f][0-9a-f])/) {
	 my $addr = hex($1);
	 my $b1 = hex($2);
	 my $b2 = hex($3);
	 my $b3 = hex($4);
	 my $b4 = hex($5);

	 while ($addr - $last_addr > 4) {
	    print "pad $1\n";
	    print OUT pack 'C*', 0, 0, 0, 0;
	    $last_addr += 4;
	    $size += 4;
	    }

	 print OUT pack 'C*', $b1, $b2, $b3, $b4;
	 $size += 4;
	 $last_addr = $addr;
	 }
      }

   close(IN);
   close(OUT);

   printf("size (hdr: %d  code: %d  total: %d)\n", 128, $size, $size+128);
   }
   



my ($infile, $outfile) = @ARGV;

create_plugin_image($infile, $outfile);
