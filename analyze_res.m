function analyze_res(res1, res2, angle1, angle2)
    disp('--- ���ս���Ա� ---');
    fprintf('�ļ�3 (%d��) -> ���: %.2f ��m, ��_p: %.1f, ��: %.1f\n', angle1, res1.thk, res1.drudeParam(1), res1.drudeParam(2));
    fprintf('�ļ�4 (%d��) -> ���: %.2f ��m, ��_p: %.1f, ��: %.1f\n', angle2, res2.thk, res2.drudeParam(1), res2.drudeParam(2));
    thk1 = res1.thk; thk2 = res2.thk;
    thk_diff_percent = abs(thk1-thk2)/mean([thk1,thk2])*100;
    fprintf('�����ǶȲ�õĺ�ȷֱ�Ϊ %.2f ��m �� %.2f ��m����Բ���Ϊ %.2f%%��\n', thk1, thk2, thk_diff_percent);
    if thk_diff_percent < 5
        fprintf('   -> ����: ���һ�������ã�����ɿ���\n');
        fprintf('   -> �Ƽ����ֵ: %.2f �� %.2f ��m\n', mean([thk1, thk2]), std([thk1, thk2]));
    else
        fprintf('   -> ����: ��Ȳ���ϴ�\n');
    end
end
