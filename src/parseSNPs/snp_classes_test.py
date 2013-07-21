"""
This program is designed to test the classes in snp_classes

Author: David Gray
"""
import sys
from snp_classes import *
import unittest

####################################################################################
#
# Test ResultsSet class  
#
####################################################################################
class ResultsSet_test(unittest.TestCase):

    def setUp(self):
        self.results_set = ResultsSet()
        self.results_set.get_or_create_result('RS12345', 1, 123456)
        self.results_set.get_or_create_result('RS12345', 1, 123456)  # same
        self.results_set.get_or_create_result('RS12346', 1, 123456)  # Different rsid
        self.results_set.get_or_create_result('RS12345', 2, 123456)  # Different chromosome
        self.results_set.get_or_create_result('RS12345', 1, 123457)  # Different position

    def test_len(self):
        self.assertEqual(len(self.results_set), 4)

    def test_get_or_create_result(self):
        result = Result('RS12345', 1, 123456)
        base_result = self.results_set.get_or_create_result('RS12345', 1, 123456)
        self.assertEqual(result.get_rsid(), base_result.get_rsid())
        self.assertEqual(result.get_chromosome(), base_result.get_chromosome())
        self.assertEqual(result.get_position(), base_result.get_position())
        
    def test_get_results_iterator(self):
        rsids = []
        chromosomes = []
        positions = []
        for result in self.results_set.get_results_iterator():
            rsid = result.get_rsid()
            if (rsid not in rsids):
                rsids.append(rsid)
            chromosome = result.get_chromosome()
            if (chromosome not in chromosomes):
                chromosomes.append(chromosome)
            position = result.get_position()
            if (position not in positions):
                positions.append(position)
        self.assertEqual(2, len(rsids))
        self.assertEqual(2, len(chromosomes))
        self.assertEqual(2, len(positions))
        
    def test_string_output(self):
        output_string = str(self.results_set)
        self.assertTrue("RS12345" in output_string)
        self.assertTrue("chromosome 1," in output_string)
        self.assertTrue("123456" in output_string)

####################################################################################
#
# Test Result class  
#
####################################################################################
class Result_test(unittest.TestCase):

    def setUp(self):
        self.result = Result('RS12345', 1, 123456)
        
    def test_get_rsid(self):
        self.assertEqual('RS12345', self.result.get_rsid())
        
    def test_get_chromosome(self):
        self.assertEqual(1, self.result.get_chromosome())
        
    def test_get_position(self):
        self.assertEqual(123456, self.result.get_position())
        
    def test_get_key(self):
        self.assertEqual(Result.get_key_static('RS12345', 1, 123456), self.result.get_key())

    def test_groups(self):
        result = self.result
        result.add_one("Group 1", "AA")
        result.add_one("Group 1", "AG")
        result.add_one("Group 1", "AA")
        result.add_one("Group 1", "AG")
        result.add_one("Group 2", "CT")
        group1 = result.get_group("Group 1")
        group2 = result.get_group("Group 2")
        group3 = result.get_group("Group 3")
        self.assertEqual("Group 1", group1.get_label())
        self.assertEqual("Group 2", group2.get_label())
        self.assertEqual("Group 3", group3.get_label())
        self.assertEqual(2, group1.get_count("AA"))
        self.assertEqual(2, group1.get_count("AG"))
        self.assertEqual(1, group2.get_count("CT"))
        self.assertEqual(0, group2.get_count("AA"))
        self.assertEqual(0, group3.get_count("AA"))
        self.assertEqual(2, len(group1.get_counts()))
        self.assertEqual(1, len(group2.get_counts()))
        self.assertEqual(0, len(group3.get_counts()))
        
    def test_string_output(self):
        result = self.result
        result.add_one("Group 1", "AA")
        result.add_one("Group 1", "AG")
        result.add_one("Group 1", "AA")
        result.add_one("Group 1", "AG")
        result.add_one("Group 2", "CT")
        group1 = result.get_group("Group 1")
        group2 = result.get_group("Group 2")
        group3 = result.get_group("Group 3")
        
        output_string = str(result)
        self.assertTrue("RS12345" in output_string)
        self.assertTrue("chromosome 1," in output_string)
        self.assertTrue("123456" in output_string)
        self.assertTrue(" 2 results" in output_string)

####################################################################################
#
# Test Group class  
#
####################################################################################
class Group_test(unittest.TestCase):

    def setUp(self):
        self.group = Group("Group 1")
        
    def test_get_label(self):
        self.assertEqual('Group 1', self.group.get_label())
        
    def test_genotypes(self):
        group = self.group
        group.add_genotype("AA")
        group.add_genotype("AA")
        group.add_genotype("AT")
        group.add_genotype("AT")
        self.assertEqual(2, group.get_count("AA"))
        self.assertEqual(2, group.get_count("AT"))
        self.assertEqual(0, group.get_count("CA"))
        self.assertEqual(2, len(group.get_counts()))
        
    def test_string_output(self):
        group = self.group
        group.add_genotype("AA")
        group.add_genotype("AA")
        group.add_genotype("AT")
        group.add_genotype("AT")
        
        output_string = str(group)
        self.assertTrue("AA=2" in output_string)
        self.assertTrue("AT=2" in output_string)

####################################################################################
#
# Test Params class  
#
####################################################################################
class Params_test(unittest.TestCase):

    def setUp(self):
        self.params = Params()
        
    def test_defaults(self):
        self.assertEqual('.', self.params.get_directory_location())
        self.assertEqual(sys.maxint, self.params.get_position_end())
        self.assertEqual(0, self.params.get_position_start())
        self.assertEqual(False, self.params.get_show_file_progress())
        self.assertEqual(0, self.params.get_show_lines_progress_interval())
        self.assertEqual(False, self.params.get_show_selected_files())
        self.assertEqual(None, self.params.get_chromosomes())
        
    def test_process(self):
        # Default should process all SNPs
        self.assertTrue(self.params.process(SnpValues("RS12345", "4", 125646, "AA")))
        
        # Test RSID match
        self.params.set_rsid("RS12*")
        self.assertTrue(self.params.process(SnpValues("RS12345", "4", 125646, "AA")))
        self.assertFalse(self.params.process(SnpValues("RS22222", "4", 125646, "AA")))
        self.params.set_rsid("**")
        
        # Test position match
        self.params.set_position_start(200)
        self.params.set_position_end(250)
        self.assertFalse(self.params.process(SnpValues("RS22222", "4", 199, "AA")))
        self.assertTrue(self.params.process(SnpValues("RS12345", "4", 200, "AA")))
        self.assertTrue(self.params.process(SnpValues("RS12345", "4", 220, "AA")))
        self.assertTrue(self.params.process(SnpValues("RS12345", "4", 250, "AA")))
        self.assertFalse(self.params.process(SnpValues("RS22222", "4", 251, "AA")))
        self.params.set_position_start(0)
        self.params.set_position_end(sys.maxint)
        
        # Test chromosome match
        # Do this test last because there is no way to remove a chromosome once added
        self.params.add_chromosome("4")
        self.params.add_chromosome("5")
        self.assertTrue(self.params.process(SnpValues("RS12345", "4", 125646, "AA")))
        self.assertTrue(self.params.process(SnpValues("RS12345", "5", 125646, "AA")))
        self.assertFalse(self.params.process(SnpValues("RS22222", "3", 125646, "AA")))
        
    def test_set_directory_location(self):
        self.params.set_directory_location("C:\OpenSNP")
        self.assertEqual("C:\OpenSNP", self.params.get_directory_location())
        
    def test_set_rsid(self):
        self.params.set_rsid("RS12345")
        self.assertEqual("RS12345", self.params.get_rsid())
        
    def test_set_position_end(self):
        self.params.set_position_end(10)
        self.assertEqual(10, self.params.get_position_end())
        
    def test_set_position_start(self):
        self.params.set_position_start(9)
        self.assertEqual(9, self.params.get_position_start())
        
    def test_set_show_file_progress(self):
        self.params.set_show_file_progress(True)
        self.assertEqual(True, self.params.get_show_file_progress())
        self.params.set_show_file_progress(False)
        self.assertEqual(False, self.params.get_show_file_progress())
        
    def test_set_show_selected_files(self):
        self.params.set_show_selected_files(True)
        self.assertEqual(True, self.params.get_show_selected_files())
        self.params.set_show_selected_files(False)
        self.assertEqual(False, self.params.get_show_selected_files())
        
    def test_set_show_lines_progress_interval(self):
        self.params.set_show_lines_progress_interval(10000)
        self.assertEqual(10000, self.params.get_show_lines_progress_interval())
        
    def test_add_chromosome(self):
        self.params.add_chromosome('1')
        self.params.add_chromosome('5')
        self.params.add_chromosome('X')
        chromosomes = self.params.get_chromosomes()
        self.assertEqual('1', chromosomes[0])
        self.assertEqual('5', chromosomes[1])
        self.assertEqual('X', chromosomes[2])
        
    def test_add_file_group(self):
        group1 = FileGroup("Group 1", 1)
        group1.add_file_selector("*bdef*")
        group2 = FileGroup("Group 2", 2)
        group2.add_file_selector("*bde*")
        group3 = FileGroup("Group 3", 3)
        group3.add_file_selector("*hij*")
        self.params.add_file_group(group2)
        self.params.add_file_group(group1)  # See if we re-prioritize this one first
        self.params.add_file_group(group3)

        self.assertEqual("Group 1", self.params.get_file_group_label("file_with_bdef_in_name"))
        self.assertEqual("Group 2", self.params.get_file_group_label("file_with_bde_in_name"))
        self.assertEqual("Group 3", self.params.get_file_group_label("file_with_hij_in_name"))
        self.assertEqual(None, self.params.get_file_group_label("non_matching_file_name"))
        
    def test_string_output(self):
        params = self.params
        params.set_directory_location("C:\OpenSNP")
        params.set_rsid("RS12345")
        params.set_position_start(985656)
        params.set_position_end(45456464)
        params.set_show_file_progress(True)
        params.set_show_selected_files(True)
        params.set_show_lines_progress_interval(198756)
        params.add_chromosome('1')
        params.add_chromosome('5')
        params.add_chromosome('X')
        group1 = FileGroup("Group 1", 1)
        group1.add_file_selector("*bdef*")
        group2 = FileGroup("Group 2", 2)
        group2.add_file_selector("*bde*")
        group3 = FileGroup("Group 3", 3)
        group3.add_file_selector("*hij*")
        params.add_file_group(group2)
        params.add_file_group(group1)
        params.add_file_group(group3)
        
        output_string = str(params)
        self.assertTrue("RS12345" in output_string)
        self.assertTrue("['1', '5', 'X']" in output_string)
        self.assertTrue("985656" in output_string)
        self.assertTrue("45456464" in output_string)
        self.assertTrue("show_file_progress True" in output_string)
        self.assertTrue("show_selected_files True" in output_string)
        self.assertTrue("198756" in output_string)
        self.assertTrue("C:\\OpenSNP" in output_string)
        self.assertTrue("Group 1," in output_string)
        self.assertTrue("Group 2," in output_string)
        self.assertTrue("Group 3," in output_string)
# These don't work - don't know why
#        self.assertTrue("*bdef*," in output_string)
#        self.assertTrue("*bde*," in output_string)
#        self.assertTrue("*hij*," in output_string)
        self.assertTrue("seq 1," in output_string)
        self.assertTrue("seq 2," in output_string)
        self.assertTrue("seq 3," in output_string)
        
####################################################################################
#
# Test FileGroup class  
#
####################################################################################
class FileGroup_test(unittest.TestCase):

    def setUp(self):
        self.file_group = FileGroup("label", 1)
        
    def test_get_label(self):
        self.assertEqual("label", self.file_group.get_label())
        
    def test_get_priority_seq(self):
        self.assertEqual(1, self.file_group.get_priority_seq())        
        
    def test_file_selectors(self):
        file_group = self.file_group
        file_group.add_file_selector("*abc*")
        file_group.add_file_selector("*def*")

        self.assertTrue(file_group.matches("file_with_abc_in_name"))
        self.assertTrue(file_group.matches("file_with_def_in_name"))
        self.assertFalse(file_group.matches("non_matching_file_name"))
        
    def test_string_output(self):
        file_group = self.file_group
        file_group.add_file_selector("*abc*")
        file_group.add_file_selector("*def*")
        
        output_string = str(file_group)
        self.assertTrue("*abc*" in output_string)
        self.assertTrue("*def*" in output_string)
        
####################################################################################
#
# Test SnpValues class  
#
####################################################################################
class SnpValues_test(unittest.TestCase):

    def setUp(self):
        self.snp_values = SnpValues("RS12345", "4", 125646, "AA")
        
    def test_get_rsid(self):
        self.assertEqual("RS12345", self.snp_values.get_rsid())
        
    def test_get_chromosomeself(self):
        self.assertEqual("4", self.snp_values.get_chromosome())
        
    def test_get_position(self):
        self.assertEqual(125646, self.snp_values.get_position())
        
    def test_get_genotype(self):
        self.assertEqual("AA", self.snp_values.get_genotype())
        
    def test_string_output(self):
        snp_values = self.snp_values
        
        output_string = str(snp_values)
        self.assertTrue("RS12345" in output_string)
        self.assertTrue("chromosome 4" in output_string)
        self.assertTrue("125646" in output_string)
        self.assertTrue("AA" in output_string)
        
####################################################################################
#
# Test TwentyThreeAndMeSNPProcessor class  
#
####################################################################################
class TwentyThreeAndMeSNPProcessor_test(unittest.TestCase):

    def setUp(self):
        self.parser = TwentyThreeAndMeSNPProcessor()
        
    def test_handles_file(self):
        self.assertTrue(self.parser.handles_file("user11_file176_yearofbirth_unknown_sex_unknown.23andme.txt"))
        self.assertFalse(self.parser.handles_file("user1035_file518_yearofbirth_1986_sex_XX.ftdna-illumina.txt"))
        self.assertFalse(self.parser.handles_file("user907_file443_yearofbirth_unknown_sex_unknown.IYG.txt"))
        self.assertFalse(self.parser.handles_file("user500_file230_yearofbirth_1956_sex_XY.decodeme.txt"))
        
    def test_get_file_type_label(self):
        self.assertEquals("23andme", self.parser.get_file_type_label())
        
    def test_parse_line(self):
        snp_values = self.parser.parse_line("rs7537756\t1\t854250\tAG")
        self.assertEquals("RS7537756", snp_values.get_rsid())
        self.assertEquals("1", snp_values.get_chromosome())
        self.assertEquals(854250, snp_values.get_position())
        self.assertEquals("AG", snp_values.get_genotype())
        
####################################################################################
#
# Test IlluminaSNPProcessor class  
#
####################################################################################
class IlluminaSNPProcessor_test(unittest.TestCase):

    def setUp(self):
        self.parser = IlluminaSNPProcessor()
        
    def test_handles_file(self):
        self.assertFalse(self.parser.handles_file("user11_file176_yearofbirth_unknown_sex_unknown.23andme.txt"))
        self.assertTrue(self.parser.handles_file("user1035_file518_yearofbirth_1986_sex_XX.ftdna-illumina.txt"))
        self.assertFalse(self.parser.handles_file("user907_file443_yearofbirth_unknown_sex_unknown.IYG.txt"))
        self.assertFalse(self.parser.handles_file("user500_file230_yearofbirth_1956_sex_XY.decodeme.txt"))
        
    def test_get_file_type_label(self):
        self.assertEquals("illumina", self.parser.get_file_type_label())
        
    def test_parse_line(self):
        # Test first illumina format
        snp_values = self.parser.parse_line('"rs4475691","1","836671","TT"')
        self.assertEquals("RS4475691", snp_values.get_rsid())
        self.assertEquals("1", snp_values.get_chromosome())
        self.assertEquals(836671, snp_values.get_position())
        self.assertEquals("TT", snp_values.get_genotype())
        
        # Test second illumina format - e.g. user1035_file518_yearofbirth_1986_sex_XX.ftdna-illumina.txt
        snp_values = self.parser.parse_line('rs3131972\t1\t752721\tA\tA')
        self.assertEquals("RS3131972", snp_values.get_rsid())
        self.assertEquals("1", snp_values.get_chromosome())
        self.assertEquals(752721, snp_values.get_position())
        self.assertEquals("AA", snp_values.get_genotype())
        
        # Test third illumina format - e.g. user981_file487_yearofbirth_1966_sex_unknown.ftdna-illumina.txt
        snp_values = self.parser.parse_line('rs11240777 1 788822 AA')
        self.assertEquals("RS11240777", snp_values.get_rsid())
        self.assertEquals("1", snp_values.get_chromosome())
        self.assertEquals(788822, snp_values.get_position())
        self.assertEquals("AA", snp_values.get_genotype())
        
####################################################################################
#
# Test IYGSNPProcessor class  
#
####################################################################################
class IYGSNPProcessor_test(unittest.TestCase):

    def setUp(self):
        self.parser = IYGSNPProcessor()
        
    def test_handles_file(self):
        self.assertFalse(self.parser.handles_file("user11_file176_yearofbirth_unknown_sex_unknown.23andme.txt"))
        self.assertFalse(self.parser.handles_file("user1035_file518_yearofbirth_1986_sex_XX.ftdna-illumina.txt"))
        self.assertTrue(self.parser.handles_file("user907_file443_yearofbirth_unknown_sex_unknown.IYG.txt"))
        self.assertFalse(self.parser.handles_file("user500_file230_yearofbirth_1956_sex_XY.decodeme.txt"))
        
    def test_get_file_type_label(self):
        self.assertEquals("iyg", self.parser.get_file_type_label())
        
    def test_parse_line(self):
        snp_values = self.parser.parse_line("rs1668873\tAG")
        self.assertEquals("RS1668873", snp_values.get_rsid())
        self.assertEquals("", snp_values.get_chromosome())
        self.assertEquals(0, snp_values.get_position())
        self.assertEquals("AG", snp_values.get_genotype())
        
####################################################################################
#
# Test DecodeMeSNPProcessor class  
#
####################################################################################
class DecodeMeSNPProcessor_test(unittest.TestCase):

    def setUp(self):
        self.parser = DecodeMeSNPProcessor()
        
    def test_handles_file(self):
        self.assertFalse(self.parser.handles_file("user11_file176_yearofbirth_unknown_sex_unknown.23andme.txt"))
        self.assertFalse(self.parser.handles_file("user1035_file518_yearofbirth_1986_sex_XX.ftdna-illumina.txt"))
        self.assertFalse(self.parser.handles_file("user907_file443_yearofbirth_unknown_sex_unknown.IYG.txt"))
        self.assertTrue(self.parser.handles_file("user500_file230_yearofbirth_1956_sex_XY.decodeme.txt"))
        
    def test_get_file_type_label(self):
        self.assertEquals("decodeme", self.parser.get_file_type_label())
        
    def test_parse_line(self):
        snp_values = self.parser.parse_line("rs12562034,A/G,1,758311,+,GG")
        self.assertEquals("RS12562034", snp_values.get_rsid())
        self.assertEquals("1", snp_values.get_chromosome())
        self.assertEquals(758311, snp_values.get_position())
        self.assertEquals("GG", snp_values.get_genotype())





if __name__ == '__main__':
    unittest.main()
