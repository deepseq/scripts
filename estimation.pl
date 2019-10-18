#!/usr/bin/perl -w
use strict;
use warnings;
use Math::Matrix;
use Getopt::Std;
use vars qw($opt_r $opt_d $opt_s $opt_w $opt_l $opt_q $opt_p $opt_o);
getopts("r:d:s:w:l:q:p:o:");

my $starttime=time;
my $runningtime=localtime(time);

my $usage="usage:\n\tperl $0 -option parameter ...\n\noptions:\n\t-r\tmatrix file\n\t-d\tmatrix from count\n\t-s\tcount matrix\n\t-w\texon file\n\t-l\tlength file\n\t-q\tread length\n\t-p\tuse matrix or not (type \"yes\" or \"no\")\n\t-o\toutput file prefix\n";
my $usagenote="Notice:\n";
my $rcfile=$opt_r || die "$usage\n$usagenote\n";
my $rcfile_s=$opt_d || die "$usage\n$usagenote\n";
my $srcfile=$opt_s || die "$usage\n$usagenote\n";
my $ewfile=$opt_w || die "$usage\n$usagenote\n";
my $elfile=$opt_l || die "$usage\n$usagenote\n";
my $readlength=$opt_q || die "$usage\n$usagenote\n";
my $pmatrix=$opt_p || die "$usage\n$usagenote\n";
my $outfile=$opt_o || die "$usage\n$usagenote\n";

open RC,"<$rcfile" or die "Cannot open infile to read.\n";
open RCS,"<$rcfile_s" or die "Cannot open infile to read.\n";
open SRC,"<$srcfile" or die "Cannot open infile to read.\n";
open EW,"<$ewfile" or die "Cannot open infile to read.\n";
open EL,"<$elfile" or die "Cannot open infile to read.\n";
open OUT,">$outfile.iso.mle.exp.txt" or die "Cannot open outfile to write.\n";

print "--$runningtime--\tData is under processing, please wait...\n\n";

my $filelinecmd="";
my @fileline=0;
my $filelinecount=0;
my $propercent=0;

sub Ltheta_iso
{
		my($ex_old,$ex_new,$prec)=@_;
		
		if(abs($ex_new-$ex_old)<$prec)
		{
				return 0;
		}
		else
		{
				return 1;
		}
}

sub resfunc_spliceread
{
		my($iso,$ex,$rl,$elm,$ewm,$rcm,$rcm_s,$srcm,$mrow,$mcolumn,$msplicemcolumn,$subp)=@_;
		my $subi=0;
		my $subj=0;
		my $subk=0;
		my @unireadnum;
		my @bothreadnum;
		my @otherreadnum;
		my $sumread=0;
		my $sumunireadnum=0;
		my $sumbothreadnum=0;
		my $sumotherreadnum=0;
		my @s_unireadnum;
		my @s_bothreadnum;
		my @s_otherreadnum;
		my $s_sumread=0;
		my $s_sumunireadnum=0;
		my $s_sumbothreadnum=0;
		my $s_sumotherreadnum=0;
		my $subres1=0;
		my $subres2=0;
		my $subres3=0;
		my $subres4=0;
		my $transsum1=0;
		my $transsum2=0;
		my $transpro1=0;
		my $transpro2=0;
		
		for($subi=0;$subi<$mrow;$subi++)
		{
				for($subj=0;$subj<$mcolumn;$subj++)
				{
						$unireadnum[$subi][$subj]=0;
						$bothreadnum[$subi][$subj]=0;
						$otherreadnum[$subi][$subj]=0;
				}
		}

		for($subj=0;$subj<$mcolumn;$subj++)
		{
				for($subi=0;$subi<$mrow;$subi++)
				{
						if($$rcm[$subi][$subj]==0)
						{
								for($subk=0;$subk<$mrow;$subk++)
								{
										if($$rcm[$subk][$subj]!=0)
										{
												$otherreadnum[$subi][$subj]=$$rcm[$subk][$subj]-$$rcm_s[$subk][$subj];
												last;
										}
								}
						}
						else
						{
								for($subk=0;$subk<$mrow;$subk++)
								{
										if($subi!=$subk)
										{
												$sumread=$sumread+$$rcm[$subk][$subj];
										}
								}
								if($sumread==0)
								{
										$unireadnum[$subi][$subj]=$$rcm[$subi][$subj]-$$rcm_s[$subi][$subj];
								}
								else
								{
										$bothreadnum[$subi][$subj]=$$rcm[$subi][$subj]-$$rcm_s[$subi][$subj];
								}
						}
						$sumread=0;
				}	
		}
		
		if($subp eq "yes")
		{
				for($subj=0;$subj<$mcolumn;$subj++)
				{
						if($$ewm[$iso][$subj]!=0)
						{
								$sumunireadnum=$sumunireadnum+$unireadnum[$iso][$subj]/$$ewm[$iso][$subj]*100;
								$sumbothreadnum=$sumbothreadnum+$bothreadnum[$iso][$subj]/$$ewm[$iso][$subj]*100;
								$sumotherreadnum=$sumotherreadnum+$otherreadnum[$iso][$subj]/$$ewm[$iso][$subj]*100;
						}
				}
		}
		else
		{
				for($subj=0;$subj<$mcolumn;$subj++)
				{
						$sumunireadnum=$sumunireadnum+$unireadnum[$iso][$subj];
						$sumbothreadnum=$sumbothreadnum+$bothreadnum[$iso][$subj];
						$sumotherreadnum=$sumotherreadnum+$otherreadnum[$iso][$subj];
				}
		}
		
		for($subi=0;$subi<$mrow;$subi++)
		{
				for($subj=0;$subj<$msplicemcolumn;$subj++)
				{
						$s_unireadnum[$subi][$subj]=0;
						$s_bothreadnum[$subi][$subj]=0;
						$s_otherreadnum[$subi][$subj]=0;
				}
		}

		for($subj=0;$subj<$msplicemcolumn;$subj++)
		{
				for($subi=0;$subi<$mrow;$subi++)
				{
						if($$srcm[$subi][$subj]==0)
						{
								for($subk=0;$subk<$mrow;$subk++)
								{
										if($$srcm[$subk][$subj]!=0)
										{
												$s_otherreadnum[$subi][$subj]=$$srcm[$subk][$subj];
												last;
										}
								}
						}
						else
						{
								for($subk=0;$subk<$mrow;$subk++)
								{
										if($subi!=$subk)
										{
												$s_sumread=$s_sumread+$$srcm[$subk][$subj];
										}
								}
								if($s_sumread==0)
								{
										$s_unireadnum[$subi][$subj]=$$srcm[$subi][$subj];
								}
								else
								{
										$s_bothreadnum[$subi][$subj]=$$srcm[$subi][$subj];
								}
						}
						$s_sumread=0;
				}	
		}
		
		for($subj=0;$subj<$msplicemcolumn;$subj++)
		{
				$s_sumunireadnum=$s_sumunireadnum+$s_unireadnum[$iso][$subj];
				$s_sumbothreadnum=$s_sumbothreadnum+$s_bothreadnum[$iso][$subj];
				$s_sumotherreadnum=$s_sumotherreadnum+$s_otherreadnum[$iso][$subj];
		}
		
		$s_sumunireadnum=$s_sumunireadnum+$sumunireadnum;
		$s_sumbothreadnum=$s_sumbothreadnum+$sumbothreadnum;
		$s_sumotherreadnum=$s_sumotherreadnum+$sumotherreadnum;
		
		for($subj=0;$subj<$mcolumn;$subj++)
		{
				$transsum1=$transsum1+$$elm[$iso][$subj];
		}
		$transpro1=1/($transsum1-$rl+1);
		
		for($subj=0;$subj<$mcolumn;$subj++)
		{
				for($subi=0;$subi<$mrow;$subi++)
				{
						if($subi!=$iso)
						{
								if($$elm[$subi][$subj]!=0)
								{
										$transsum2=$transsum2+$$elm[$subi][$subj];
										last;
								}
						}
				}
		}
		$transpro2=1/($transsum2-$rl+1);
		
		$subres1=(2*$s_sumunireadnum+$s_sumotherreadnum+$s_sumbothreadnum)*$transpro2-($s_sumunireadnum+$s_sumbothreadnum)*$transpro1;
		$subres2=4*$s_sumunireadnum*$transpro2*(($s_sumunireadnum+$s_sumotherreadnum+$s_sumbothreadnum)*($transpro1-$transpro2));
		$subres3=((2*$s_sumunireadnum+$s_sumotherreadnum+$s_sumbothreadnum)*$transpro2-($s_sumunireadnum+$s_sumbothreadnum)*$transpro1)**2;
		$subres4=2*(($s_sumunireadnum+$s_sumotherreadnum+$s_sumbothreadnum)*($transpro2-$transpro1));
		
		if($subres4 != 0)
		{
				return ($subres1-sqrt($subres2+$subres3))/$subres4;
		}
		else
		{
				return -1;
		}
}

sub resfind_spliceread
{
		my($iso,$ex,$rl,$elm,$ewm,$rcm,$rcm_s,$srcm,$mrow,$mcolumn,$msplicemcolumn,$subp)=@_;
		my $subi=0;
		my $subj=0;
		my $subk=0;
		my @unireadnum;
		my @bothreadnum;
		my @otherreadnum;
		my $sumread=0;
		my $sumunireadnum=0;
		my $sumbothreadnum=0;
		my $sumotherreadnum=0;
		my @s_unireadnum;
		my @s_bothreadnum;
		my @s_otherreadnum;
		my $s_sumread=0;
		my $s_sumunireadnum=0;
		my $s_sumbothreadnum=0;
		my $s_sumotherreadnum=0;
		my $subres1=0;
		my $subres2=0;
		my $subres3=0;
		my $subres4=0;
		my $transsum1=0;
		my $transsum2=0;
		my $transpro1=0;
		my $transpro2=0;
		
		for($subi=0;$subi<$mrow;$subi++)
		{
				for($subj=0;$subj<$mcolumn;$subj++)
				{
						$unireadnum[$subi][$subj]=0;
						$bothreadnum[$subi][$subj]=0;
						$otherreadnum[$subi][$subj]=0;
				}
		}

		for($subj=0;$subj<$mcolumn;$subj++)
		{
				for($subi=0;$subi<$mrow;$subi++)
				{
						if($$rcm[$subi][$subj]==0)
						{
								for($subk=0;$subk<$mrow;$subk++)
								{
										if($$rcm[$subk][$subj]!=0)
										{
												$otherreadnum[$subi][$subj]=$$rcm[$subk][$subj]-$$rcm_s[$subk][$subj];
												last;
										}
								}
						}
						else
						{
								for($subk=0;$subk<$mrow;$subk++)
								{
										if($subi!=$subk)
										{
												$sumread=$sumread+$$rcm[$subk][$subj];
										}
								}
								if($sumread==0)
								{
										$unireadnum[$subi][$subj]=$$rcm[$subi][$subj]-$$rcm_s[$subi][$subj];
								}
								else
								{
										$bothreadnum[$subi][$subj]=$$rcm[$subi][$subj]-$$rcm_s[$subi][$subj];
								}
						}
						$sumread=0;
				}	
		}
		
		if($subp eq "yes")
		{
				for($subj=0;$subj<$mcolumn;$subj++)
				{
						if($$ewm[$iso][$subj]!=0)
						{
								$sumunireadnum=$sumunireadnum+$unireadnum[$iso][$subj]/$$ewm[$iso][$subj]*100;
								$sumbothreadnum=$sumbothreadnum+$bothreadnum[$iso][$subj]/$$ewm[$iso][$subj]*100;
								$sumotherreadnum=$sumotherreadnum+$otherreadnum[$iso][$subj]/$$ewm[$iso][$subj]*100;
						}
				}
		}
		else
		{
				for($subj=0;$subj<$mcolumn;$subj++)
				{
						$sumunireadnum=$sumunireadnum+$unireadnum[$iso][$subj];
						$sumbothreadnum=$sumbothreadnum+$bothreadnum[$iso][$subj];
						$sumotherreadnum=$sumotherreadnum+$otherreadnum[$iso][$subj];
				}
		}
		
		for($subi=0;$subi<$mrow;$subi++)
		{
				for($subj=0;$subj<$msplicemcolumn;$subj++)
				{
						$s_unireadnum[$subi][$subj]=0;
						$s_bothreadnum[$subi][$subj]=0;
						$s_otherreadnum[$subi][$subj]=0;
				}
		}

		for($subj=0;$subj<$msplicemcolumn;$subj++)
		{
				for($subi=0;$subi<$mrow;$subi++)
				{
						if($$srcm[$subi][$subj]==0)
						{
								for($subk=0;$subk<$mrow;$subk++)
								{
										if($$srcm[$subk][$subj]!=0)
										{
												$s_otherreadnum[$subi][$subj]=$$srcm[$subk][$subj];
												last;
										}
								}
						}
						else
						{
								for($subk=0;$subk<$mrow;$subk++)
								{
										if($subi!=$subk)
										{
												$s_sumread=$s_sumread+$$srcm[$subk][$subj];
										}
								}
								if($s_sumread==0)
								{
										$s_unireadnum[$subi][$subj]=$$srcm[$subi][$subj];
								}
								else
								{
										$s_bothreadnum[$subi][$subj]=$$srcm[$subi][$subj];
								}
						}
						$s_sumread=0;
				}	
		}
		
		for($subj=0;$subj<$msplicemcolumn;$subj++)
		{
				$s_sumunireadnum=$s_sumunireadnum+$s_unireadnum[$iso][$subj];
				$s_sumbothreadnum=$s_sumbothreadnum+$s_bothreadnum[$iso][$subj];
				$s_sumotherreadnum=$s_sumotherreadnum+$s_otherreadnum[$iso][$subj];
		}
		
		$s_sumunireadnum=$s_sumunireadnum+$sumunireadnum;
		$s_sumbothreadnum=$s_sumbothreadnum+$sumbothreadnum;
		$s_sumotherreadnum=$s_sumotherreadnum+$sumotherreadnum;
		
		$s_sumunireadnum=$s_sumunireadnum+$s_sumbothreadnum*$$ex[$iso];
		$s_sumotherreadnum=$s_sumotherreadnum+$s_sumbothreadnum*(1-$$ex[$iso]);
		$s_sumbothreadnum=0;
		
		for($subj=0;$subj<$mcolumn;$subj++)
		{
				$transsum1=$transsum1+$$elm[$iso][$subj];
		}
		$transpro1=1/($transsum1-$rl+1);
		
		for($subj=0;$subj<$mcolumn;$subj++)
		{
				for($subi=0;$subi<$mrow;$subi++)
				{
						if($subi!=$iso)
						{
								if($$elm[$subi][$subj]!=0)
								{
										$transsum2=$transsum2+$$elm[$subi][$subj];
										last;
								}
						}
				}
		}
		$transpro2=1/($transsum2-$rl+1);
		
		$subres1=(2*$s_sumunireadnum+$s_sumotherreadnum+$s_sumbothreadnum)*$transpro2-($s_sumunireadnum+$s_sumbothreadnum)*$transpro1;
		$subres2=4*$s_sumunireadnum*$transpro2*(($s_sumunireadnum+$s_sumotherreadnum+$s_sumbothreadnum)*($transpro1-$transpro2));
		$subres3=((2*$s_sumunireadnum+$s_sumotherreadnum+$s_sumbothreadnum)*$transpro2-($s_sumunireadnum+$s_sumbothreadnum)*$transpro1)**2;
		$subres4=2*(($s_sumunireadnum+$s_sumotherreadnum+$s_sumbothreadnum)*($transpro2-$transpro1));
		
		if($subres4 != 0)
		{
				return ($subres1-sqrt($subres2+$subres3))/$subres4;
		}
		else
		{
				return -1;
		}
}

sub mle_lrt
{
		my($iso,$ex,$rl,$elm,$ewm,$rcm,$rcm_s,$srcm,$mrow,$mcolumn,$msplicemcolumn,$subp)=@_;
		my $subi=0;
		my $subj=0;
		my $subk=0;
		my @unireadnum;
		my @bothreadnum;
		my @otherreadnum;
		my $sumread=0;
		my $sumunireadnum=0;
		my $sumbothreadnum=0;
		my $sumotherreadnum=0;
		my @s_unireadnum;
		my @s_bothreadnum;
		my @s_otherreadnum;
		my $s_sumread=0;
		my $s_sumunireadnum=0;
		my $s_sumbothreadnum=0;
		my $s_sumotherreadnum=0;
		my $subres1=0;
		my $subres2=0;
		my $subres3=0;
		my $subres4=0;
		my $transsum1=0;
		my $transsum2=0;
		my $transpro1=0;
		my $transpro2=0;
		my $theta_random=0;
		my $theta_uni=0;
		my $theta_other=0;
		my $theta_both=0;
		my $lrt_value=0;
		my $p_num=0;
		
		for($subi=0;$subi<$mrow;$subi++)
		{
				for($subj=0;$subj<$mcolumn;$subj++)
				{
						$unireadnum[$subi][$subj]=0;
						$bothreadnum[$subi][$subj]=0;
						$otherreadnum[$subi][$subj]=0;
				}
		}

		for($subj=0;$subj<$mcolumn;$subj++)
		{
				for($subi=0;$subi<$mrow;$subi++)
				{
						if($$rcm[$subi][$subj]==0)
						{
								for($subk=0;$subk<$mrow;$subk++)
								{
										if($$rcm[$subk][$subj]!=0)
										{
												$otherreadnum[$subi][$subj]=$$rcm[$subk][$subj]-$$rcm_s[$subk][$subj];
												last;
										}
								}
						}
						else
						{
								for($subk=0;$subk<$mrow;$subk++)
								{
										if($subi!=$subk)
										{
												$sumread=$sumread+$$rcm[$subk][$subj];
										}
								}
								if($sumread==0)
								{
										$unireadnum[$subi][$subj]=$$rcm[$subi][$subj]-$$rcm_s[$subi][$subj];
								}
								else
								{
										$bothreadnum[$subi][$subj]=$$rcm[$subi][$subj]-$$rcm_s[$subi][$subj];
								}
						}
						$sumread=0;
				}	
		}
		
		if($subp eq "yes")
		{
				for($subj=0;$subj<$mcolumn;$subj++)
				{
						if($$ewm[$iso][$subj]!=0)
						{
								$sumunireadnum=$sumunireadnum+$unireadnum[$iso][$subj]/$$ewm[$iso][$subj]*100;
								$sumbothreadnum=$sumbothreadnum+$bothreadnum[$iso][$subj]/$$ewm[$iso][$subj]*100;
								$sumotherreadnum=$sumotherreadnum+$otherreadnum[$iso][$subj]/$$ewm[$iso][$subj]*100;
						}
				}
		}
		else
		{
				for($subj=0;$subj<$mcolumn;$subj++)
				{
						$sumunireadnum=$sumunireadnum+$unireadnum[$iso][$subj];
						$sumbothreadnum=$sumbothreadnum+$bothreadnum[$iso][$subj];
						$sumotherreadnum=$sumotherreadnum+$otherreadnum[$iso][$subj];
				}
		}
		
		for($subi=0;$subi<$mrow;$subi++)
		{
				for($subj=0;$subj<$msplicemcolumn;$subj++)
				{
						$s_unireadnum[$subi][$subj]=0;
						$s_bothreadnum[$subi][$subj]=0;
						$s_otherreadnum[$subi][$subj]=0;
				}
		}

		for($subj=0;$subj<$msplicemcolumn;$subj++)
		{
				for($subi=0;$subi<$mrow;$subi++)
				{
						if($$srcm[$subi][$subj]==0)
						{
								for($subk=0;$subk<$mrow;$subk++)
								{
										if($$srcm[$subk][$subj]!=0)
										{
												$s_otherreadnum[$subi][$subj]=$$srcm[$subk][$subj];
												last;
										}
								}
						}
						else
						{
								for($subk=0;$subk<$mrow;$subk++)
								{
										if($subi!=$subk)
										{
												$s_sumread=$s_sumread+$$srcm[$subk][$subj];
										}
								}
								if($s_sumread==0)
								{
										$s_unireadnum[$subi][$subj]=$$srcm[$subi][$subj];
								}
								else
								{
										$s_bothreadnum[$subi][$subj]=$$srcm[$subi][$subj];
								}
						}
						$s_sumread=0;
				}	
		}
		
		for($subj=0;$subj<$msplicemcolumn;$subj++)
		{
				$s_sumunireadnum=$s_sumunireadnum+$s_unireadnum[$iso][$subj];
				$s_sumbothreadnum=$s_sumbothreadnum+$s_bothreadnum[$iso][$subj];
				$s_sumotherreadnum=$s_sumotherreadnum+$s_otherreadnum[$iso][$subj];
		}
		
		$s_sumunireadnum=$s_sumunireadnum+$sumunireadnum;
		$s_sumbothreadnum=$s_sumbothreadnum+$sumbothreadnum;
		$s_sumotherreadnum=$s_sumotherreadnum+$sumotherreadnum;
		
		$theta_uni=$s_sumunireadnum;
		$theta_other=$s_sumotherreadnum;
		$theta_both=$s_sumbothreadnum;
		
		$s_sumunireadnum=$s_sumunireadnum+$s_sumbothreadnum*$$ex[$iso];
		$s_sumotherreadnum=$s_sumotherreadnum+$s_sumbothreadnum*(1-$$ex[$iso]);
		$s_sumbothreadnum=0;
			
		for($subj=0;$subj<$mcolumn;$subj++)
		{
				$transsum1=$transsum1+$$elm[$iso][$subj];
		}
		$transpro1=1/($transsum1-$rl+1);
		
		for($subj=0;$subj<$mcolumn;$subj++)
		{
				for($subi=0;$subi<$mrow;$subi++)
				{
						if($subi!=$iso)
						{
								if($$elm[$subi][$subj]!=0)
								{
										$transsum2=$transsum2+$$elm[$subi][$subj];
										last;
								}
						}
				}
		}
		$transpro2=1/($transsum2-$rl+1);
		
		if(($$ex[$iso]>0)&&($$ex[$iso]<1))
		{
				for($subi=0;$subi<1000;$subi++)
				{
						while($theta_random==0)
						{
								$theta_random=rand(1);					
						}
						
						$theta_uni=$theta_uni+$theta_both*$theta_random;
						$theta_other=$theta_other+$theta_both*(1-$theta_random);
						
						$lrt_value=2*($theta_uni*log($transpro1*$theta_random)+$theta_other*log($transpro2*(1-$theta_random))-$s_sumunireadnum*log($transpro1*$$ex[$iso])-$s_sumotherreadnum*log($transpro2*(1-$$ex[$iso])));
						
						if($lrt_value<3.84)
						{
								$p_num++;
						}
						
						$theta_uni=$theta_uni-$theta_both*$theta_random;
						$theta_other=$theta_other-$theta_both*(1-$theta_random);	
						
						$theta_random=0;					
				}
		}
		else
		{
				$p_num=0;
		}

		return $p_num/1000;
}

my @ReadCountMatrix=<RC>;
my @ReadCountMatrixSplice=<RCS>;
my @SpliceReadCountMatrix=<SRC>;
my @ExonWeightMatrix=<EW>;
my @ExonLengthMatrix=<EL>;

my @tmpstrRC="";
my @tmpstrRCS="";
my @tmpstrSRC="";
my @tmpstrEW="";
my @tmpstrEL="";
my $genename="";
my @accnum="";

my $i=0;
my $j=0;
my $k=0;
my $rows=0;
my $columns=0;
my $splicecolumns=0;
my $m=0;
my $n=0;
my $l=0;

my @readcount;
my @readcountsplice;
my @splicereadcount;
my @exonweight;
my @exonlength;
my @s_theta_old;
my @s_theta_new;
my @s_p_value;

my $donext=1;

for($i=0;$i<($#ReadCountMatrix+1);$i++)
{
		chomp $ReadCountMatrix[$i];
		chomp $ReadCountMatrixSplice[$i];
		chomp $SpliceReadCountMatrix[$i];
		chomp $ExonWeightMatrix[$i];
		chomp $ExonLengthMatrix[$i];
		@tmpstrRC=split("\t",$ReadCountMatrix[$i]);
		@tmpstrRCS=split("\t",$ReadCountMatrixSplice[$i]);
		@tmpstrSRC=split("\t",$SpliceReadCountMatrix[$i]);
		@tmpstrEW=split("\t",$ExonWeightMatrix[$i]);
		@tmpstrEL=split("\t",$ExonLengthMatrix[$i]);
		if($genename eq "")
		{
				$genename=$tmpstrRC[0];
				$accnum[$l]=$tmpstrRC[1];
				$l++;
				for($j=2;$j<($#tmpstrRC+1);$j++)
				{
						$readcount[$k][$j-2]=$tmpstrRC[$j];
						$readcountsplice[$k][$j-2]=$tmpstrRCS[$j];
						$exonweight[$k][$j-2]=$tmpstrEW[$j];
						$exonlength[$k][$j-2]=$tmpstrEL[$j];
				}
				for($j=2;$j<($#tmpstrSRC+1);$j++)
				{
						$splicereadcount[$k][$j-2]=$tmpstrSRC[$j];
				}
				$k++;
				$columns=$#tmpstrRC-1;
				$splicecolumns=$#tmpstrSRC-1;
		}
		else
		{
				if($genename eq $tmpstrRC[0])
				{
						for($j=2;$j<($#tmpstrRC+1);$j++)
						{
								$readcount[$k][$j-2]=$tmpstrRC[$j];
								$readcountsplice[$k][$j-2]=$tmpstrRCS[$j];
								$exonweight[$k][$j-2]=$tmpstrEW[$j];
								$exonlength[$k][$j-2]=$tmpstrEL[$j];
						}
						for($j=2;$j<($#tmpstrSRC+1);$j++)
						{
								$splicereadcount[$k][$j-2]=$tmpstrSRC[$j];
						}
						$k++;
						
						$accnum[$l]=$tmpstrRC[1];
						$l++;
						
						if($i==$#ReadCountMatrix)
						{
								$rows=$k;
						
								for($m=0;$m<$rows;$m++)
								{
										$s_theta_new[$m]=0.1;
										$s_theta_old[$m]=0.2;
								}
								
								if($rows==1)
								{
										$s_theta_new[0]=1;
										$s_p_value[0]=1;
								}
								else
								{		
										for($m=0;$m<$rows;$m++)
										{												
												$s_theta_new[$m]=resfunc_spliceread($m,\@s_theta_old,$readlength,\@exonlength,\@exonweight,\@readcount,\@readcountsplice,\@splicereadcount,$rows,$columns,$splicecolumns,$pmatrix);
												while(Ltheta_iso($s_theta_old[$m],$s_theta_new[$m],0.001))
												{
														$s_theta_old[$m]=$s_theta_new[$m];
														$s_theta_new[$m]=resfind_spliceread($m,\@s_theta_old,$readlength,\@exonlength,\@exonweight,\@readcount,\@readcountsplice,\@splicereadcount,$rows,$columns,$splicecolumns,$pmatrix);
												}
												$s_p_value[$m]=mle_lrt($m,\@s_theta_new,$readlength,\@exonlength,\@exonweight,\@readcount,\@readcountsplice,\@splicereadcount,$rows,$columns,$splicecolumns,$pmatrix);
										}
								}
						
								for($m=0;$m<$rows;$m++)
								{
										print OUT "$genename\t$accnum[$m]\t$s_p_value[$m]\t$s_theta_new[$m]\n";
								}
										
								for($m=0;$m<$rows;$m++)
								{
										for($n=0;$n<$columns;$n++)
										{
												$readcount[$m][$n]="";
												$readcountsplice[$m][$n]="";
												$exonweight[$m][$n]="";
												$exonlength[$m][$n]="";
										}
								}
								for($m=0;$m<$rows;$m++)
								{
										for($n=0;$n<$splicecolumns;$n++)
										{
												$splicereadcount[$m][$n]="";
										}
								}
								
								@s_theta_old="";
								@s_theta_new="";
								@s_p_value="";

								$genename="";
								@accnum="";
								$l=0;
								$k=0;
								$rows=0;
								$columns=0;
								$splicecolumns=0;						
						}
				}
				else
				{
						$rows=$k;
						
						for($m=0;$m<$rows;$m++)
						{
								$s_theta_new[$m]=0.1;
								$s_theta_old[$m]=0.2;
						}
						
						if($rows==1)
						{
								$s_theta_new[0]=1;
								$s_p_value[0]=1;
						}		
						else
						{
								for($m=0;$m<$rows;$m++)
								{									
										$s_theta_new[$m]=resfunc_spliceread($m,\@s_theta_old,$readlength,\@exonlength,\@exonweight,\@readcount,\@readcountsplice,\@splicereadcount,$rows,$columns,$splicecolumns,$pmatrix);
										while(Ltheta_iso($s_theta_old[$m],$s_theta_new[$m],0.001))
										{
												$s_theta_old[$m]=$s_theta_new[$m];
												$s_theta_new[$m]=resfind_spliceread($m,\@s_theta_old,$readlength,\@exonlength,\@exonweight,\@readcount,\@readcountsplice,\@splicereadcount,$rows,$columns,$splicecolumns,$pmatrix);
										}
										$s_p_value[$m]=mle_lrt($m,\@s_theta_new,$readlength,\@exonlength,\@exonweight,\@readcount,\@readcountsplice,\@splicereadcount,$rows,$columns,$splicecolumns,$pmatrix);										
								}
						}
						
						for($m=0;$m<$rows;$m++)
						{
								print OUT "$genename\t$accnum[$m]\t$s_p_value[$m]\t$s_theta_new[$m]\n";
						}
						
						for($m=0;$m<$rows;$m++)
						{
								for($n=0;$n<$columns;$n++)
								{
										$readcount[$m][$n]="";
										$readcountsplice[$m][$n]="";
										$exonweight[$m][$n]="";
										$exonlength[$m][$n]="";
								}
						}
						for($m=0;$m<$rows;$m++)
						{
								for($n=0;$n<$splicecolumns;$n++)
								{
										$splicereadcount[$m][$n]="";
								}
						}
						
						@s_theta_old="";
						@s_theta_new="";
						@s_p_value="";
						
						$genename="";
						@accnum="";
						$l=0;
						$k=0;
						$rows=0;
						$columns=0;			
						$splicecolumns=0;			
						
						$genename=$tmpstrRC[0];
						$accnum[$l]=$tmpstrRC[1];
						$l++;
						for($j=2;$j<($#tmpstrRC+1);$j++)
						{
								$readcount[$k][$j-2]=$tmpstrRC[$j];
								$readcountsplice[$k][$j-2]=$tmpstrRCS[$j];
								$exonweight[$k][$j-2]=$tmpstrEW[$j];
								$exonlength[$k][$j-2]=$tmpstrEL[$j];
						}
						for($j=2;$j<($#tmpstrSRC+1);$j++)
						{
								$splicereadcount[$k][$j-2]=$tmpstrSRC[$j];
						}
						$k++;
						
						$columns=$#tmpstrRC-1;
						$splicecolumns=$#tmpstrSRC-1;
						
						if($i==$#ReadCountMatrix)
						{
								$rows=$k;
						
								for($m=0;$m<$rows;$m++)
								{
										$s_theta_new[$m]=0.1;
										$s_theta_old[$m]=0.2;
								}
									
								if($rows==1)
								{
										$s_theta_new[0]=1;
										$s_p_value[0]=1;
								}
								else
								{		
										for($m=0;$m<$rows;$m++)
										{													
												$s_theta_new[$m]=resfunc_spliceread($m,\@s_theta_old,$readlength,\@exonlength,\@exonweight,\@readcount,\@readcountsplice,\@splicereadcount,$rows,$columns,$splicecolumns,$pmatrix);
												while(Ltheta_iso($s_theta_old[$m],$s_theta_new[$m],0.001))
												{
														$s_theta_old[$m]=$s_theta_new[$m];
														$s_theta_new[$m]=resfind_spliceread($m,\@s_theta_old,$readlength,\@exonlength,\@exonweight,\@readcount,\@readcountsplice,\@splicereadcount,$rows,$columns,$splicecolumns,$pmatrix);
												}
												$s_p_value[$m]=mle_lrt($m,\@s_theta_new,$readlength,\@exonlength,\@exonweight,\@readcount,\@readcountsplice,\@splicereadcount,$rows,$columns,$splicecolumns,$pmatrix);
										}
								}
						
								for($m=0;$m<$rows;$m++)
								{
										print OUT "$genename\t$accnum[$m]\t$s_p_value[$m]\t$s_theta_new[$m]\n";
								}
										
								for($m=0;$m<$rows;$m++)
								{
										for($n=0;$n<$columns;$n++)
										{
												$readcount[$m][$n]="";
												$readcountsplice[$m][$n]="";
												$exonweight[$m][$n]="";
												$exonlength[$m][$n]="";
										}
								}
								for($m=0;$m<$rows;$m++)
								{
										for($n=0;$n<$splicecolumns;$n++)
										{
												$splicereadcount[$m][$n]="";
										}
								}
								
								@s_theta_old="";
								@s_theta_new="";
								@s_p_value="";
								
								$genename="";
								@accnum="";
								$l=0;
								$k=0;
								$rows=0;
								$columns=0;	
								$splicecolumns=0;					
						}						
				}
		}
		
		@tmpstrRC="";
		@tmpstrSRC="";
		@tmpstrEW="";
		@tmpstrEL="";
		
		if($filelinecount==int(($#ReadCountMatrix+1)/100))
		{
				$runningtime=localtime(time);
				$propercent++;
				$filelinecount=0;
				print STDERR "--$runningtime--\tRead Processing Percentage:$propercent%\r";
		}
		else
		{
				$filelinecount++;
		}
}

my $endtime=time;
my $costtime=$endtime-$starttime;
print "Now calculation is finished.\n";
print "running time: $costtime\n";

close(OUT);
close(EL);
close(EW);
close(SRC);
close(RC);
